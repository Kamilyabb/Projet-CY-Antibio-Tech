import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Créer les dossiers nécessaires
os.makedirs("output", exist_ok=True)
os.makedirs("images", exist_ok=True)

# Étape 1 : Lire les données avec le bon séparateur
data = pd.read_csv("input/data_real.csv", sep=';')

# Étape 2 : Afficher les colonnes disponibles
print("Colonnes disponibles :", data.columns)

# Étape 3 : Normaliser les noms de colonnes
data.columns = data.columns.str.strip().str.replace(' ', '_').str.lower()
print("Colonnes après normalisation :", data.columns)

# Étape 4 : Afficher les valeurs uniques de certaines colonnes importantes
print("Valeurs uniques dans 'sample_type' :", data['sample_type'].unique())
print("Valeurs uniques dans 'treatment' :", data['treatment'].unique())

# Étape 5 : Filtrage des données fécales
if 'sample_type' in data.columns and 'mouse_id' in data.columns and 'counts_live_bacteria_per_wet_g' in data.columns:
    fecal_data = data[data['sample_type'] == 'fecal'][['mouse_id', 'treatment', 'experimental_day', 'counts_live_bacteria_per_wet_g']]
    print("Données fécales après filtrage :")
    print(fecal_data.head())
    fecal_data.to_csv("output/fecal_data.csv", index=False)
    print("Données fécales filtrées avec succès.")
else:
    print("Erreur : Certaines colonnes nécessaires sont introuvables dans les données.")

# Étape 6 : Création des courbes multi-segments (fécal)
if not fecal_data.empty:
    print("Nombre de lignes dans fecal_data :", len(fecal_data))
    plt.figure(figsize=(10, 6))
    for mouse, group in fecal_data.groupby('mouse_id'):
        print(f"Traitement de la souris {mouse}...")
        treatment = group['treatment'].iloc[0]
        color = 'blue' if treatment == 'placebo' else 'orange'
        plt.plot(group['experimental_day'], group['counts_live_bacteria_per_wet_g'], 
                 label=f"{treatment} - Mouse {mouse}" if mouse == 1 else "", 
                 alpha=0.6, color=color)

    plt.title("Fecal Live Bacteria")
    plt.xlabel("Washout day")
    plt.ylabel("log10(live bacteria/wet g)")
    plt.legend(['Placebo', 'ABX'])
    plt.grid(True)
    plt.savefig("images/fecal_graph.png")
    print("Graphique fécal sauvegardé : 'images/fecal_graph.png'")
else:
    print("Pas de données fécales disponibles pour créer le graphique.")

# Étape 7 : Filtrage des données pour le graphique violon (cécal)
if 'sample_type' in data.columns and 'experimental_day' in data.columns and 'counts_live_bacteria_per_wet_g' in data.columns:
    cecal_data = data[(data['sample_type'] == 'cecal') & (data['experimental_day'] == 0)][['treatment', 'counts_live_bacteria_per_wet_g']]
    print("Données cécales après filtrage :")
    print(cecal_data.head())
    cecal_data.to_csv("output/cecal_data.csv", index=False)
    print("Données cécales filtrées avec succès.")
else:
    print("Erreur : Données cécales non disponibles.")

# Étape 8 : Création du graphique violon (cécal)
if not cecal_data.empty:
    plt.figure(figsize=(8, 6))
    sns.violinplot(data=cecal_data, x='treatment', y='counts_live_bacteria_per_wet_g', 
                   palette={'placebo': 'blue', 'ABX': 'orange'})
    plt.title("Cecal Live Bacteria")
    plt.xlabel("Treatment")
    plt.ylabel("log10(live bacteria/wet g)")
    plt.savefig("images/cecal_graph.png")
    print("Graphique cécal sauvegardé : 'images/cecal_graph.png'")
else:
    print("Pas de données cécales disponibles pour créer le graphique.")

# Étape 9 : Filtrage des données pour le graphique violon (iléon)
if 'sample_type' in data.columns and 'experimental_day' in data.columns and 'counts_live_bacteria_per_wet_g' in data.columns:
    ileal_data = data[(data['sample_type'] == 'ileal') & (data['experimental_day'] == 0)][['treatment', 'counts_live_bacteria_per_wet_g']]
    print("Données iléales après filtrage :")
    print(ileal_data.head())
    ileal_data.to_csv("output/ileal_data.csv", index=False)
    print("Données iléales filtrées avec succès.")
else:
    print("Erreur : Données iléales non disponibles.")

# Étape 10 : Création du graphique violon (iléon)
if not ileal_data.empty:
    plt.figure(figsize=(8, 6))
    sns.violinplot(data=ileal_data, x='treatment', y='counts_live_bacteria_per_wet_g', 
                   palette={'placebo': 'blue', 'ABX': 'orange'})
    plt.title("Ileal Live Bacteria")
    plt.xlabel("Treatment")
    plt.ylabel("log10(live bacteria/wet g)")
    plt.savefig("images/ileal_graph.png")
    print("Graphique iléal sauvegardé : 'images/ileal_graph.png'")
else:
    print("Pas de données iléales disponibles pour créer le graphique.")

print("Tous les fichiers et graphiques ont été générés avec succès.")
