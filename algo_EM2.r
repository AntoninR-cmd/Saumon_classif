# Importation des données
dta <- read.table("mixture.csv", header = FALSE, col.names = "X", dec = ",")
summary(dta)
head(dta)

# Affichage de l'histogramme
hist(dta$X, breaks = 30, col = "lightblue", main = "Histogramme des données", xlab = "Valeurs")

# Les paramètres initiaux sont choisis en fonction de l'histogramme
# Cependant, dans la boucle ils sont calculés  à partir d'une attribution initiale des groupes
# Selon l'histogramme on différencie les deux groupes avec une séparation à la valeur 25

dta$tau_1 <- ifelse(dta$X < 25, 1, 0)
dta$tau_2 <- 1-dta$tau_1

# On fixe un seuil arbitraire de 10^⁻3
seuil <- 10^-3

# On défnit différence à 1 pour lancer la boucle et itération à 0 pour compter le nombre d'itérations
diff <- 1
iteration <- 0

# on définit l'ancienne vraisemblance à -Inf pour la première itération
acn_vrais <- -Inf

# On boucle tant que la différence de vraisemblance entre deux itérations est supérieure au seuil
while (diff > seuil) {
    iteration <- iteration + 1

    # Etape M on calcule les paramètres du modèle selon les groupes obtenus
    n1 <- sum(dta$tau_1)
    n2 <- sum(dta$tau_2)

    # On calcule les proportions de chaque groupe
    pi_1 <- n1/nrow(dta)
    pi_2 <- 1 - pi_1

    # On actualise les nouveaux paramètres du modèle selon les groupes obtenus
    mu1 <- sum(dta$tau_1*dta$X)/n1
    sd1 <- sqrt(sum(dta$tau_1*(dta$X - mu1)^2)/n1)

    mu2 <- sum(dta$tau_2*dta$X)/n2
    sd2 <- sqrt(sum(dta$tau_2*(dta$X - mu2)^2)/n2)

    # On calcule la proba selon la formule de Bayes
    # tau = proba du groupe k * vraisemblance du x_i dans groupe k / somme des vraisemblances de tous les groupes
    # On calcule la densité de chaque observation dans chaque groupe avec dnorm
    # Etape E
    dta$dens_1 <- dnorm(dta$X, mean = mu1, sd = sd1)
    dta$dens_2 <- dnorm(dta$X, mean = mu2, sd = sd2)

    # On calcule les numérateurs de la formule de Bayes pour chaque observation et chaque groupe
    dta$num_1 <- (pi_1 * dta$dens_1)
    dta$num_2 <- (pi_2 * dta$dens_2)

    # On calcule de dénominateur de la formule de Bayes pour chaque observation soit la somme des proba
    dta$somme <- dta$num_1 + dta$num_2

    # On calcule les tau pour chaque observation et chaque groupe selon la fromule de Bayes
    dta$tau_1 <- dta$num_1 / dta$somme
    dta$tau_2 <- dta$num_2 / dta$somme

    # On calcule la vraisemblance somme des log de somme des proba
    nvl_vrais <- sum(log(dta$somme))

    # On calcule la vraisemblance du modèle selon les paramètres actuels
    diff <- abs(nvl_vrais - acn_vrais)

    # Si la diff est > à 10-3, on actualise l'ancienne vraisemblance pour la prochaine itération
    if (diff > seuil) {
        acn_vrais <- nvl_vrais
    }
}

# Après la convergence des tau on affecte aux groupes selon la plus grande proba
dta$groupe <- ifelse(dta$tau_1 > dta$tau_2, 1, 2)

print(paste("Nombre d'itérations :", iteration))
print(paste("Groupe 1 : mu =", round(mu1, 2), ", sd =", round(sd1, 2), ", proportion =", round(pi_1, 2)))
print(paste("Groupe 2 : mu =", round(mu2, 2), ", sd =", round(sd2, 2), ", proportion =", round(pi_2, 2)))
print(table(dta$groupe))
