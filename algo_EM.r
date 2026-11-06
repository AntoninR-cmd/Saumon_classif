# Algorithme EM codé à la main
# Mélange de 2 lois gaussiennes : on cherche mu1, sd1, mu2, sd2 et la
# proportion p du groupe 1 (le groupe 2 a la proportion 1 - p).


# ---------------------------------------------------------------------------
# 1. Importation des données
# ---------------------------------------------------------------------------

# Le fichier n'a pas d'en-tête et utilise la virgule comme séparateur décimal
# (32,90563243 = 32.9), d'où read.csv2.
mixture <- read.csv2("mixture.csv", header = FALSE)
head(mixture)

x <- mixture$V1   # le vecteur à étudier
n <- length(x)    # nombre d'individus
n


# ---------------------------------------------------------------------------
# 2. Histogramme
# ---------------------------------------------------------------------------

hist(x, breaks = 30, freq = FALSE, col = "lightgrey",
     main = "Histogramme du vecteur", xlab = "x")

# On voit deux bosses : une étroite vers 21 et une plus large vers 33.


# ---------------------------------------------------------------------------
# 3. Initialisation des paramètres à partir de l'histogramme
# ---------------------------------------------------------------------------

# Groupe 1 : pic vers 21, presque tout entre 17 et 25 -> écart-type d'environ 2
mu1 <- 21
sd1 <- 2

# Groupe 2 : pic vers 33, presque tout entre 25 et 41 -> écart-type d'environ 4
mu2 <- 33
sd2 <- 4

# On ne sait pas encore quel groupe est le plus gros : on part sur 50/50
p <- 0.5


# ---------------------------------------------------------------------------
# Densité d'une loi normale écrite à la main (pas de dnorm)
# ---------------------------------------------------------------------------

densite_normale <- function(x, mu, sd) {
  1 / (sd * sqrt(2 * pi)) * exp(-(x - mu)^2 / (2 * sd^2))
}


# ---------------------------------------------------------------------------
# 4. à 6. Boucle EM
# ---------------------------------------------------------------------------

seuil <- 0.001          # arrêt quand la vraisemblance bouge de moins de 0.001
ancienne_vrais <- -Inf  # pas de vraisemblance avant la première itération
difference <- Inf       # pour entrer dans la boucle
iteration <- 0
historique <- c()       # on garde la vraisemblance de chaque itération

while (difference > seuil) {

  iteration <- iteration + 1

  # --- 4. Calcul de la vraisemblance -------------------------------------

  # Densité de chaque individu s'il venait du groupe 1
  f1 <- densite_normale(x, mu1, sd1)

  # Densité de chaque individu s'il venait du groupe 2
  f2 <- densite_normale(x, mu2, sd2)

  # Densité du mélange pour chaque individu : chaque groupe est pondéré
  # par sa proportion
  f_melange <- p * f1 + (1 - p) * f2

  # Les individus sont indépendants : la vraisemblance est le produit
  # des densités. Un produit de 650 nombres petits donne 0 pour l'ordinateur,
  # donc on passe au log, qui transforme le produit en somme.
  vrais <- sum(log(f_melange))

  # On garde la valeur pour tracer la convergence à la fin
  historique <- c(historique, vrais)

  # Écart avec l'itération précédente (EM fait toujours monter la vraisemblance)
  difference <- vrais - ancienne_vrais

  # La valeur actuelle devient l'ancienne pour le tour suivant
  ancienne_vrais <- vrais

  # --- 5. Attribution des individus (étape E) -----------------------------

  # Probabilité que chaque individu appartienne au groupe 1, sachant sa valeur
  # (formule de Bayes). Chaque individu est partagé entre les deux groupes :
  # il compte pour t1 dans le groupe 1 et pour t2 dans le groupe 2.
  t1 <- p * f1 / f_melange
  t2 <- 1 - t1

  # --- 6. Calcul des nouveaux paramètres (étape M) ------------------------

  # Proportion du groupe 1 : la part moyenne des individus qui lui revient
  p <- sum(t1) / n

  # Moyennes : moyennes de x pondérées par les probabilités d'appartenance
  mu1 <- sum(t1 * x) / sum(t1)
  mu2 <- sum(t2 * x) / sum(t2)

  # Écarts-types : même chose avec les écarts au carré à la nouvelle moyenne
  sd1 <- sqrt(sum(t1 * (x - mu1)^2) / sum(t1))
  sd2 <- sqrt(sum(t2 * (x - mu2)^2) / sum(t2))

  cat("Itération", iteration,
      "| log-vraisemblance =", round(vrais, 4),
      "| mu1 =", round(mu1, 2), "sd1 =", round(sd1, 2),
      "| mu2 =", round(mu2, 2), "sd2 =", round(sd2, 2),
      "| p =", round(p, 3), "\n")
}


# ---------------------------------------------------------------------------
# Résultats
# ---------------------------------------------------------------------------

cat("\nConvergence atteinte en", iteration, "itérations\n")
resultats <- data.frame(groupe     = c(1, 2),
                        moyenne    = c(mu1, mu2),
                        ecart_type = c(sd1, sd2),
                        proportion = c(p, 1 - p))
resultats

# Évolution de la log-vraisemblance au fil des itérations
plot(historique, type = "b", pch = 19,
     xlab = "Itération", ylab = "Log-vraisemblance",
     main = "Convergence de l'algorithme EM")

# Histogramme avec les deux gaussiennes estimées et le mélange
hist(x, breaks = 30, freq = FALSE, col = "lightgrey",
     main = "Mélange estimé par EM", xlab = "x")
grille <- seq(min(x), max(x), length.out = 500)
lines(grille, p * densite_normale(grille, mu1, sd1), col = "blue", lwd = 2)
lines(grille, (1 - p) * densite_normale(grille, mu2, sd2), col = "red", lwd = 2)
lines(grille, p * densite_normale(grille, mu1, sd1) +
        (1 - p) * densite_normale(grille, mu2, sd2), lwd = 2, lty = 2)
legend("topright", legend = c("Groupe 1", "Groupe 2", "Mélange"),
       col = c("blue", "red", "black"), lty = c(1, 1, 2), lwd = 2)

# Classe finale : chaque individu va dans le groupe où sa probabilité
# d'appartenance est la plus forte
classe <- ifelse(t1 > 0.5, 1, 2)
table(classe)
