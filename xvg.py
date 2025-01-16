import numpy as np

# .xvg dosyasını okuma
def read_xvg(file_path):
    time = []
    coulomb_energy = []
    lj_energy = []
    
    with open(file_path, 'r') as file:
        for line in file:
            # Yorum satırlarını atla
            if line.startswith('@') or line.startswith('#'):
                continue
            # Verileri ayır
            columns = line.split()
            if len(columns) >= 3:
                time.append(float(columns[0]))
                coulomb_energy.append(float(columns[1]))
                lj_energy.append(float(columns[2]))
    
    return np.array(time), np.array(coulomb_energy), np.array(lj_energy)

# Enerji bileşenlerinin toplamını ve hata paylarını hesaplama
def calculate_total_energy(coulomb_energy, lj_energy):
    # Toplam enerji
    total_energy = coulomb_energy + lj_energy
    # Enerji bileşenlerinin standart sapmalarını hesapla
    std_coul = np.std(coulomb_energy)
    std_lj = np.std(lj_energy)
    # Toplam hata payı
    total_error = np.sqrt(std_coul**2 + std_lj**2)
    
    # Ortalama toplam enerji
    mean_total_energy = np.mean(total_energy)
    
    return mean_total_energy, total_error

# Dosya yolunu belirleyin
file_path = 'interaction_energy.xvg'

# Verileri oku
time, coulomb_energy, lj_energy = read_xvg(file_path)

# Toplam enerjiyi hesapla
mean_total_energy, total_error = calculate_total_energy(coulomb_energy, lj_energy)

# Sonuçları yazdır
print(f"Toplam Enerji: {mean_total_energy:.2f} ± {total_error:.2f} kJ/mol")

