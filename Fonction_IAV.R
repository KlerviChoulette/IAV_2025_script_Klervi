################################################################################
########################### FONCTION POUR L'IAV ################################


remove_zero = function(PopN){
  
  #Quelles lignes ont des zéro à remplacer
  IndexZero = which(PopN == 0)
  
  #Si la moyenne de la population est 0 alors on mets une petite valeur
  if (mean(PopN) == 0) {
    OffsetVal = 1e-17
  } else {
    #Indice des lignes sans zéro
    IndexNonZero = which(PopN != 0)
    # 1% de la moyenne des popvalues sans les zéros
    OffsetVal = mean(PopN[IndexNonZero]) * 0.01
  }
  
  #Ajouter cette valeur aux lignes qui comportent un zéro
  if (length(IndexZero) > 0) {
    PopN = PopN + OffsetVal
  }
  return(PopN)
}



fit_gam_model = function(SmoothParm, PopNLog, YearPop){
  
  # Calcul du gam sur le log des populations
  s = mgcv::s
  model = mgcv::gam(PopNLog ~ s(YearPop,
                                 k = SmoothParm),
                    #method = "REML"
                    # +s(YearPop,
                    #     k = SmoothParm, bs = "re")
                    , fx = TRUE)
  
  
  ### Validation du gam ###
  
  rsd = residuals(model)
  
  # Calcul du gam sur les résidus
  s = mgcv::s
  modelres = mgcv::gam(rsd ~ s(YearPop, 
                                k = length(PopN), bs = "cs"), gamma = 1.4)

  ##Vérifier l'autocorréaltion
  #acf(residuals(model))
  
  ##Autres types de validation du GAM
  # plot(model,
  #      pages = 1, residuals = TRUE, all.terms = TRUE, shade = TRUE,
  #      shade.col = 2)
  #readline(prompt = "Press any key to continue")
  #mgcv::gam.check(model)
  #readline(prompt = "Press any key to continue")

  if ((abs(sum(modelres$edf) - 1)) < 0.01) {
    
    #Si le GAM est validé alors on fait le prédict des nouvelles valeurs de popvalue
    
    PopNInt = predict(model, data.frame(YearPop = YearPopInt), se.fit = T)
    #Puis exponentiel pour revenir au vrai valeur
    PopNInt_se = exp(PopNInt$se.fit)
    PopNInt = exp(PopNInt$fit)
    
    return(list("PopNInt" = PopNInt,"PopNInt_se" = PopNInt_se, "Modele" = model))
  }
}


fit_gamm_model = function(SmoothParm, PopNLog, YearPop){
  
  # Calcul du gam sur le log des populations
  s = mgcv::s
  model = mgcv::gam(PopNLog ~ s(YearPop,
                                k = SmoothParm)
                    +s(YearPop,
                        k = SmoothParm, bs = "re")
                    , fx = TRUE)
  
  
  ### Validation du gam ###
  
  rsd = residuals(model)
  
  # Calcul du gam sur les résidus
  s = mgcv::s
  modelres = mgcv::gam(rsd ~ s(YearPop, 
                               k = length(PopN), bs = "cs"), gamma = 1.4)
  
  ##Vérifier l'autocorréaltion
  #acf(residuals(model))
  
  ##Autres types de validation du GAM
  # plot(model,
  #      pages = 1, residuals = TRUE, all.terms = TRUE, shade = TRUE,
  #      shade.col = 2)
  #readline(prompt = "Press any key to continue")
  #mgcv::gam.check(model)
  #readline(prompt = "Press any key to continue")
  
  if ((abs(sum(modelres$edf) - 1)) < 0.01) {
    
    #Si le GAM est validé alors on fait le prédict des nouvelles valeurs de popvalue
    
    PopNInt = predict(model, data.frame(YearPop = YearPopInt), se.fit = T)
    #Puis exponentiel pour revenir au vrai valeur
    PopNInt_se = exp(PopNInt$se.fit)
    PopNInt = exp(PopNInt$fit)
    
    return(list("PopNInt" = PopNInt,"PopNInt_se" = PopNInt_se))
  }
}


fit_linear_model = function(YearPop, PopNLog){
  model = lm(PopNLog ~ YearPop)
  PopNInt = predict(model, data.frame(YearPop = YearPopInt), se.fit = T)
  PopNInt_se = exp(PopNInt$se.fit)
  PopNInt = exp(PopNInt$fit)
  return(list("PopNInt" = PopNInt, "PopNInt_se" = PopNInt_se))
}


fit_chain_model = function(YearPopInt, YearPop,PopN){
  
  # Apply the default approach (Chain)
  PopNInt = matrix(NA, 1, length(YearPopInt))
  PopNInt_se = matrix(NA, 1, length(YearPopInt))
  for (K in 1:length(YearPopInt)) {
    k = which(YearPop == YearPopInt[K])
    if (length(k) > 0) {
      PopNInt[K] = PopN[k]
    } else {
      # find the previous value
      YearStart = YearPopInt[K]
      YearStart = YearStart - 1
      k = which(YearPop == YearStart)
      while (length(k) == 0) {
        YearStart = YearStart - 1
        k = which(YearPop == YearStart)
      }
      PopNStart = PopN[k]
      # find the next value
      YearEnd = YearPopInt[K]
      YearEnd = YearEnd + 1
      k = which(YearPop == YearEnd)
      while (length(k) == 0) {
        YearEnd = YearEnd + 1
        k = which(YearPop == YearEnd)
      }
      PopNEnd = PopN[k]
      # Calculate the interpolated value
      PopNInt[K] = PopNStart * ((PopNEnd / PopNStart)^((YearPopInt[K] -
                                                           YearStart) / (YearEnd - YearStart)))
    }
  }
  return(list("PopNInt" = PopNInt, "PopNInt_se" = PopNInt_se))
}



calcul_lambda_pop = function(PopLambda, JIndex, FinalYear,
                              PopN, InitialYear, REF_YEAR) {
  # Année de début des calculs :
    #l'année qui suit la première année d'observation (car on ne peut pas calculer lambda sans t-1)
  StartYear = min(YearPopInt) + 1
  
  # Initialisation :
    #On fixe 1 pour la première valeur de lambda (pas de croissance possible pour la première année)
  #JIndex
  PopLambda[JIndex, 1] = 1  #JIndex = l'indice de la population
  
  #Boucle sur les années pour calculer lambda
  
  for (K in StartYear:FinalYear) {
    
    # On vérifie qu'on a bien des données pour K et pour l'année précédente (K-1)
    if (!is.na(PopN[K - InitialYear + 1]) && !is.na(PopN[K - InitialYear])) {
      # Calcul du taux de croissance
      PopLambda[JIndex, K - REF_YEAR + 1] = PopN[K - InitialYear + 1] - PopN[K - InitialYear]
    } else {
      # Si une des valeurs est manquante → NA
      PopLambda[JIndex, K - REF_YEAR + 1] = NA
    }
  }
  return(PopLambda)
}


calcul_lambda_species = function(InitialYear, FinalYear, PopLambda,
                                 SpeciesIndex, LAMBDA_MAX, LAMBDA_MIN, REF_YEAR){
  EndYear = FinalYear - REF_YEAR + 1
  
  YearIndex = 1
  PopLambda_k = NULL
  
  #pour récuperer le nombre de population utiliser pour chaque lambda d'espèces
  nb_pop_year_lambda = rep(NA, EndYear)
  if(length(PopLambda) > 0){
  while (YearIndex < EndYear +1 ) {
    
    # Pour chaque année, s'il y a au moins une valeur (non NA)
    PopLambda_k = which(!is.na(PopLambda[, YearIndex]))
    nb_pop_year_lambda[YearIndex] = length(PopLambda_k) # =- AJOUT
    
    if (length(PopLambda_k) > 0) {
      
      # Only get those lambdas that are not 'NA'
      PopLambdaTemp = PopLambda[PopLambda_k, YearIndex]
      
      # Fine which of these are less than our max
      IndexTemp = which(PopLambdaTemp <= LAMBDA_MAX) # -2
      IndexTempBad_max = which(PopLambdaTemp > LAMBDA_MAX) # 2

      # If we have some...
      if (length(IndexTemp) > 0) {
        
        # Extract them as PopLambdaTemp1
        PopLambdaTemp1 = PopLambdaTemp[IndexTemp] # -2
        
        # Get those values that are also more then our min
        
        IndexTemp = which(PopLambdaTemp1 >= LAMBDA_MIN) # no

        IndexTempBad_min = which(PopLambdaTemp1 < LAMBDA_MIN) # yes, -2
        
        # If there are some...
        if (length(IndexTemp) > 0) {
          # Then set species lambda to be their average...
          
          ### ON CHERCHE A DIMINUER LES EXTREMES AVEC LAMBDA_MAX et MIN. 
          #### Donc si CAP_LAMBDAS = TRUE alors on fait la moyenne si les 
          #### valeurs correctes et sur les extrêmes remplacés par LAMBDA_MAX et MIN
          if (CAP_LAMBDAS) {
            SpeciesLambda[SpeciesIndex, YearIndex] = mean(c(PopLambdaTemp1[IndexTemp],
                                                  rep(LAMBDA_MAX, length(IndexTempBad_max)),
                                                  rep(LAMBDA_MIN, length(IndexTempBad_min))))
          } else {
            SpeciesLambda[SpeciesIndex, YearIndex] = mean(PopLambdaTemp1[IndexTemp])
          }
        } else {
          # If there are no good lambas more than min *or* less than max, we may still have to combine them
          # Otherwise, if we have lambdas less than our max, but not more then min, set sp. av. to be NA
          if (CAP_LAMBDAS) {
            # SpeciesLambda[SpeciesIndexI, YearIndex] = LAMBDA_MIN
            SpeciesLambda[SpeciesIndex, YearIndex] = mean(c(rep(LAMBDA_MAX, length(IndexTempBad_max)), rep(LAMBDA_MIN, length(IndexTempBad_min))))
          } else {
            SpeciesLambda[SpeciesIndex, YearIndex] = NA
          }
        }
      } else {
        # Otherwise, if we have no lambdas less than our max, set sp. av. to be NA
        if (CAP_LAMBDAS) {
          SpeciesLambda[SpeciesIndex, YearIndex] = LAMBDA_MAX
        } else {
          SpeciesLambda[SpeciesIndex, YearIndex] = NA
        }
      }
    } else {
      # If all values are 'NA' then set species lambda to be NA
      SpeciesLambda[SpeciesIndex, YearIndex] = NA
    }
    YearIndex = YearIndex + 1
  }
  }
  return(list("SpeciesLambda" = SpeciesLambda, "nb_pop_year_lambda" = nb_pop_year_lambda))
}



calcul_dtemp = function(SpeciesLambda, CAP_LAMBDAS){
  DTemp = matrix(0, 1, dim(SpeciesLambda)[2])
  # For each year
  for (YearIndex in 1:dim(SpeciesLambda)[2]) {
    # Get data for this year 'YearIndex'
    YearData = SpeciesLambda[, YearIndex]
    
    # Find populations that have data
    if (!CAP_LAMBDAS) {
      Index = which(YearData != -1)
    } else {
      Index = which(!is.na(YearData))
    }
    
    # If there are some populations
    if (length(Index) > 0) {
      # DTemp is mean lambda for those populations
      DTemp[YearIndex] = mean(YearData[Index])
      nb_sp_year_dtemp[YearIndex] = length(Index)

    } else {
      DTemp[YearIndex] = NA
    }
  }
  #colnames(DTemp) = InitialYear:FinalYear
  return(list("DTemp" = DTemp, "nb_sp_year_dtemp" = nb_sp_year_dtemp))
}

calcul_index = function(DTemp, DSize, REF_YEAR, InitialYear){
  # Calculate LPI
  # Calcul du nombre de NA à ajouter au début
  
  nb_NA = InitialYear - REF_YEAR
  if (nb_NA != 0 & ncol(DTemp) < DSize) {
    DTemp = t(as.matrix(c(rep(NA, nb_NA), DTemp)))
  }
  
  # Initialiser l'indice
  I = matrix(0, DSize)
  I[1] = 1
  
  # For each year indexed as 'J'
  for (J in 2:DSize) {
    # Initialisation
    DT = 0
    DI = 0
    
    # Si la valeur pour l'année J est valide, on l'utilise
    if (!is.na(DTemp[J]) && DTemp[J] != -99) {
      DT = DTemp[J]
      DI = 1
    } else {
      warning(sprintf("Year %d is NA or flagged as -99\n", J))
    }
    
    # Calcul du nouvel indice
    if (DI > 0) {
      if (is.na(I[J - 1]) || I[J - 1] == -99) {
        I[J] = 1 * 10^(DT / DI)
        warning(sprintf("**** [Year %d] Previous year data missing, assuming '1' **** \n", J))
      } else {
        I[J] = I[J - 1] * 10^(DT / DI)
      }
    } else {
      I[J] = I[J - 1] * 10^(0)
      warning(sprintf("Year %d has no valid lambda. Assuming lambda of 0\n", J))
    }
  }
  return(I)
}

calcul_index_2 = function(DTemp, DSize, REF_YEAR, InitialYear){
  # Calculate LPI
  # Calcul du nombre de NA à ajouter au début
  
  nb_NA = InitialYear - REF_YEAR
  if (nb_NA != 0 & ncol(DTemp) < DSize) {
    DTemp = t(as.matrix(c(rep(NA, nb_NA), DTemp)))
  }
  
  # Initialiser l'indice
  I = matrix(0, DSize)
  I[1] = 1
  
  # For each year indexed as 'J'
  for (J in 2:DSize) {
    # Initialisation
    DT = 0
    DI = 0
    
    # Si la valeur pour l'année J est valide, on l'utilise
    if (!is.na(DTemp[J]) && DTemp[J] != -99) {
      DT = DTemp[J]
      DI = 1
    } else {
      warning(sprintf("Year %d is NA or flagged as -99\n", J))
    }
    
    # Calcul du nouvel indice
    if (DI > 0) {
      if (is.na(I[J - 1]) || I[J - 1] == -99) {
        I[J] = 1 * 10^(DT / DI)
        warning(sprintf("**** [Year %d] Previous year data missing, assuming '1' **** \n", J))
      } else {
        I[J] = I[J - 1] * 10^(DT / DI)
      }
    } else {
      I[J] = NA
      warning(sprintf("Year %d has no valid lambda. Assuming lambda of 0\n", J))
    }
  }
  return(I)
}

calcul_bootstrap_LPI = function(SpeciesLambda, DSize, CAP_LAMBDAS, 
                                InitialYear, REF_YEAR){
  
  nb_NA = InitialYear - REF_YEAR
  
  # Initialiser l'indice
  BootI = matrix(0, DSize)
  BootI[1] = 1
  BootDTemp = matrix(0, DSize)
  BootDTemp[1] = 1
  
  for (J in 2:DSize) {
    # Récupérer les valeurs de lambda à l'année J
    SpeciesLambda_2 = SpeciesLambda[1:dim(SpeciesLambda)[1], J]
    
    # Nettoyer les valeurs manquantes
    SpeciesLambdaVal = na.omit(SpeciesLambda_2)
    
    # Bootstrap : échantillonnage avec remise
    BootVal = sample(SpeciesLambdaVal, replace = TRUE)
    
    # Appliquer filtre selon CAP_LAMBDAS
    if (!CAP_LAMBDAS) {
      Index = which(BootVal != -1)
    } else {
      Index = which(!is.na(BootVal))
    }
    
    if (length(Index) > 0) {
      # Moyenne des valeurs bootstrapées valides
      DT = mean(BootVal[Index])
      DI = 1  # Un seul groupe
    } else {
      DT = 0
      DI = 0
    }
    
    # Mise à jour de l’indice
    if (DI == 0) {
      BootI[J] = NA
    } else {
      if (is.na(BootI[J - 1])) {
        BootI[J] = 1 * 10^(DT)
        BootDTemp[J] = DT
      } else {
        BootI[J] = BootI[J - 1] * 10^(DT)
        BootDTemp[J] = DT
      }
    }
  }
  
  return(list("BootI" = BootI, "BootDTemp" = BootDTemp))
}

calcul_IC = function(BOOT_STRAP_SIZE, DSize, SpeciesLambda, 
                     InitialYear, REF_YEAR){
  
  # Create matrix for bootstrap indices
  BootI = matrix(0, BOOT_STRAP_SIZE, DSize)
  BootDTemp = matrix(0, BOOT_STRAP_SIZE, DSize)
  BootIFlag = matrix(0, 1, BOOT_STRAP_SIZE)
  
  # # Converted to parallel *********
  # BootI = foreach::foreach(Loop = 1:BOOT_STRAP_SIZE) %op% {
  #   calcul_bootstrap_lpi(SpeciesLambdaArray, DSize, CAP_LAMBDAS)
  # 
  # }
  
  for (Loop in 1:BOOT_STRAP_SIZE) {
    BootI_vect = calcul_bootstrap_LPI(SpeciesLambda, DSize, CAP_LAMBDAS, 
                                       InitialYear, REF_YEAR)
    BootDTemp_vect = BootI_vect$BootDTemp
    BootI_vect = BootI_vect$BootI
      
    
    # Ajouter BootI_vect comme ligne à la matrice BootI (en utilisant l'indice Loop)
    BootI[Loop, ] = BootI_vect
    BootDTemp[Loop, ] =BootDTemp_vect
  }
  
  # # Combine list of vectors from foreach into matrix
  # BootI = do.call(cbind, BootI)
  # # Transpose matrix as each bootstrap loop is a column and we'd like them to be rows
  # BootI = t(BootI)
  
  
  CIx = matrix(0, DSize, 2)
  CIx[1, 1] = 1
  CIx[1, 2] = 1
  CI_DTemp = matrix(0, DSize, 2)
  CI_DTemp[1, 1] = 1
  CI_DTemp[1, 2] = 1
  # Estimate confidence intervals using the bootstapped indicies
  for (J in 2:DSize) {
    # If this is a valid index year for this group
    if (valid_index_years[J] & all(!is.na(BootI[, J]))) {
      # Get the data
      BootIVal = BootI[, J]
      BootDTempVal =  BootDTemp[, J]
      # RF: this was used in original, now bootstrap only samples from valid data
      # Index = which(BootIFlag != 1)
      # BootIVal = BootI[Index, J]
      
      CIx[J, 1] = quantile(BootIVal, 0.025, names = FALSE)
      CIx[J, 2] = quantile(BootIVal, 0.975, names = FALSE)
      CI_DTemp[J, 1] = quantile(BootDTempVal, 0.025, names = FALSE)
      CI_DTemp[J, 2] = quantile(BootDTempVal, 0.975, names = FALSE)
    } else {
      # If we don't have an index for this year, we shouldn't have
      CIx[J, 1] = 1
      CIx[J, 2] = 1
      CI_DTemp[J, 1] = 1
      CI_DTemp[J, 2] = 1
    }
    
  }
  return(list("IC" = CIx, "Bootstrap" = BootI, "BootDTemp"= BootDTemp, "CI_DTemp" = CI_DTemp))
}



## Pour récuperer le taux de croissance du LPI
calcul_lambdas_LPI <- function(lpi) {
  N = length(lpi)
  lambdas = matrix(0, N)
  lambdas[1] = 1
  for (k in 2:N) {
    lambdas[k] = log10(lpi[k]) - log10(lpi[k-1])
  }
  return(lambdas)
}





# calculate the species leverage
calcul_leverage_LPI = function(SpeciesLambdaArray, DSize){
  
  I <- numeric(DSize)
  I[1] <- 1  # L'indice de base est 1
  
  SpeciesLambda_lev = matrix(ncol = DSize) 
  # Pour chaque année (à partir de la 2e)
  for (J in 2:DSize) {
    
    # Récupère les lambdas pour l'année J
    SpeciesLambdaVal <- SpeciesLambdaArray[, J]
    
    # Supprime les valeurs manquantes
    SpeciesLambdaVal <- SpeciesLambdaVal[!is.na(SpeciesLambdaVal)]
    
    # Supprime les valeurs flaguées comme invalides (-1)
    Index <- which(SpeciesLambdaVal != -1)
    #SpeciesLambda_lev[, J] = SpeciesLambdaVal[Index]
    
    # S'il reste des valeurs valides
    if (length(Index) > 0) {
      mean_lambda <- mean(SpeciesLambdaVal[Index])
      
      # Calcul de l'indice LPI (variation géométrique en base 10)
      if (is.na(I[J - 1])) {
        I[J] <- 1 * 10^mean_lambda
      } else {
        I[J] <- I[J - 1] * 10^mean_lambda
      }
    } else {
      # Pas de données valides → valeur manquante
      I[J] <- NA
    }
  }
  
  return(list("Indice" = I,"Valide_speciesLambda" = SpeciesLambdaVal))
  
}




calcul_indicateur <- function(data, var_groupe, data_SpeciesLambda_2002, CAP_LAMBDAS,
                              FinalYear, InitialYear, sNames, data_all) {
  
  asso_SPNames_Number <- tibble(
    binomial = sNames,
    number = seq_along(sNames)
  )
  
  res_Ifinal <- list()
  res_leverage <- list()
  nb_sp_max <- list()
  
  for (val in unique(data[[var_groupe]])) {
    cat("==> ", var_groupe, " : ", val, "\n")
    
    sp_group <- data %>%
      filter(.data[[var_groupe]] == val) %>%
      distinct(binomial) %>%
      pull(binomial)
    
    num_sp <- asso_SPNames_Number %>%
      filter(binomial %in% sp_group) %>%
      pull(number)
    
    if (length(num_sp) <= 1) next
    
    # Calcul principal
    res_temp <- calcul_dtemp(data_SpeciesLambda_2002[num_sp, ], CAP_LAMBDAS)
    DTemp <- res_temp$DTemp
    nb_sp_max[[val]] <- max(res_temp$nb_sp_year_dtemp)
    
    PLOT_MAX <- FinalYear
    REF_YEAR <- 2002
    DSize <- PLOT_MAX - REF_YEAR + 1
    
    Ifinal <- calcul_index(DTemp, DSize, REF_YEAR, InitialYear) %>%
      as_tibble() %>%
      mutate(LPI_final = V1, year = REF_YEAR:PLOT_MAX) %>%
      select(-V1)
    
    #Intervalle de confiance
    BOOT_STRAP_SIZE <- 1000
    IC <- calcul_IC(BOOT_STRAP_SIZE, DSize, data_SpeciesLambda_2002[num_sp, ],
                    InitialYear, REF_YEAR)$IC %>%
      as_tibble() %>%
      rename(IC_low = V1, IC_high = V2) %>%
      mutate(year = REF_YEAR:PLOT_MAX)
    
    Ifinal <- left_join(Ifinal, IC, by = "year")
    res_Ifinal[[val]] <- Ifinal
    
    #Leverage
    overall_lambdas <- calcul_lambdas_LPI(Ifinal$LPI_final)
    if(dim(data_SpeciesLambda_2002[num_sp, ])[1] <3){
      leverage_results <- lapply(seq_along(num_sp), function(i) {
        DTemp_lev <- calcul_dtemp(matrix(data_SpeciesLambda_2002[num_sp[-i], ], nrow = length(num_sp)-1), CAP_LAMBDAS)$DTemp
        lev_I <- calcul_index(DTemp_lev, DSize, 2002, 2002)
        list(
          lev_I = lev_I,
          diff = calcul_lambdas_LPI(lev_I) - overall_lambdas,
          sp = sp_group[i]
        )
      })
    }else{
      leverage_results <- lapply(seq_along(num_sp), function(i) {
        DTemp_lev <- calcul_dtemp(data_SpeciesLambda_2002[num_sp[-i], ], CAP_LAMBDAS)$DTemp
        lev_I <- calcul_index(DTemp_lev, DSize, 2002, 2002)
        list(
          lev_I = lev_I,
          diff = calcul_lambdas_LPI(lev_I) - overall_lambdas,
          sp = sp_group[i]
        )
      })}
    
    
    leverage_df <- bind_rows(lapply(leverage_results, function(l) {
      tibble(
        binomial = l$sp,
        year = seq(REF_YEAR, REF_YEAR + DSize - 1),
        LPI = as.numeric(l$lev_I)
      )
    }))
    
    res_leverage[[val]] <- leverage_df %>%
      left_join(
        data_all %>% select(binomial, all_of(var_groupe)) %>% distinct(binomial, .keep_all = TRUE),
        by = "binomial"
      )
  }
  
  list(
    Ifinal = imap_dfr(res_Ifinal, ~ mutate(.x, !!var_groupe := .y)),
    leverage = imap_dfr(res_leverage, ~ mutate(.x, !!var_groupe := .y)),
    nb_sp_max = nb_sp_max
  )
}
