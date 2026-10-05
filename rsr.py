import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats

# 1. Rastgele üretici sabitleme
np.random.seed(101)

# 2. Sentetik fits verisi ve header oluşturma
image_size = (200, 200)
# Arkaplan gürültüsü (Mean = 100 ADU, Std= 5 ADU)
background = np.random.normal(loc=100.0, scale=10.0, size=image_size)

# Görüş alanona birkaç yıldız ekleyelim
image_data = background.copy()
star_centers = [(50,60), (120, 140), (80, 150)] # (y, x) pikselleri
star_fluxes = [1500, 3000, 800] # Yıldızların tepe akıları

y_grid, x_grid = np.indices(image_size)
for (cy, cx), flux in zip(star_centers, star_fluxes):
    # Yıldız profili (PSF - Point Spread Function)
    star = flux * np.exp(-((x_grid - cx)**2 + (y_grid - cy)**2) / (2* 2.0**2))
    image_data += star

# WCS (Dünya Koordinat Sistemi) Bilgilerini Header'a yazma
wcs = WCS(naxis=2)
wcs.wcs.crpix = [100, 100]  # Referans piksel (x, y)
wcs.wcs.cdelt = [-0.0001, 0.0001] # Piksel boyutu (derece/piksel)
wcs.wcs.crval = [280.375, -0.091] # Galaktik merkeze yakın RA, Derece
wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"] #kOORDİNAT TİPİ

header = wcs.to_header()
header['EXPTIME'] = 30.0 # Poz Süresi ( saniye)
header['FILTER'] = 'V' # V Bandı filtresi

# FITS dosyasını diske kaydetme
hdu = fits.PrimaryHDU(data=image_data, header = header)
hdu.writeto('synthetic_star_field.fits', overwrite=True)
print("-> FITS dosyası başarıyla oluşturuldu: synthetic_star_field.fits")

# 3. ADIM: VERİ ANALİZİ VE İŞLEME PIPELINE'I

# A. FITS Verisini Okuma
with fits.open('synthetic_star_field.fits') as hdul:
    fits_data = hdul[0].data
    fits_header = hdul[0].header
    fits_wcs = WCS(fits_header)

# B. Arka Plan Tahmini (Sigma Clipping)
mean, median, std = sigma_clipped_stats(fits_data, sigma=3.0)
data_subtracted = fits_data - median
print(f"-> Arka Plan Seviyesi (Median): {median:.2f}, ADU, Gürültü (Std): {std:.2f} ADU")

# C. Apertür Fotometrisi (İkinci Yıldız İçin: cx=140, cy=120)
target_x, target_y = 140, 120
r_aperture = 5 # Apertür yarıçapı (piksel)

# Apertür içindeki pikselleri maskeleme
dist_from_center =np.sqrt((x_grid - target_x)**2 + (y_grid - target_y)**2)
aperture_mask = dist_from_center <= r_aperture

#Toplam net akı hesabı (Arka plan çıkarılmış veriden)
total_flux = np.sum(data_subtracted[aperture_mask])

# D. WCS Dönüşümü (Piksel -> Gökyüzü Koordinatları)
ra, dec = fits_wcs.pixel_to_world_values(target_x, target_y)

# E. Kadir hesabı ve Toz sönümlemesi Temizliği
# Temsili 0 Noktası
zp = 25.0 # Zero Point
m_obs= zp - 2.5 * np.log10(total_flux)

# Yıldızın doğrultusundaki renk fazlalığı E(B-V) =0.4 ise
E_BV = 0.4
Rv = 3.1
A_V = Rv * E_BV
m_intrinsic = m_obs - A_V # Sönümsüz (gerçek) kadir

print("\n" + "="*40)
print(f" HEDEF YILDIZ ANALİZ SONUÇLARI")
print("="*40)
print(f"Piksel Konumu: X={target_x}, Y={target_y}")
print(f"Gökyüzü Koordinatı: RA={ra:.5f} deg, DEC={dec:.5f} deg")
print(f"Ölçülen Net Akı: {total_flux:.1f} ADU")
print(f"Gözlenen Kadir (m_V): {m_obs:.3f} mag")
print(f"Toz Sönümlemesi (A_V): {A_V:.2f} mag (E(B-V)={E_BV})")
print(f"Düzeltilmiş Kadir    : {m_intrinsic:.3f} mag")
print("="*40)

# 4. Görselleştirme
plt.figure(figsize=(9, 7))
# Logaritmil ölçeklendirme ile FITS görüntüsü çizdirme
vmin, vmax = np.percentile(fits_data, [1,98])
plt.subplot(projection = fits_wcs)
plt.imshow(fits_data, origin='lower', cmap='viridis', vmin=vmin, vmax=vmax)
plt.colorbar(label='Akı (ADU)')

# Apertür halkası çizme (Hedef yıldız üzerine)
circle = plt.Circle((target_x, target_y), r_aperture, color='red', fill=False, linewidth=2, label='Fotometrik Apertür')
plt.gca().add_patch(circle)

plt.xlabel('Sağ Açıklık (RA)')
plt.ylabel('Dik Açıklık (DEC)')
plt.title('FITS VERİSİ ÜZERİNDE WCS VE APERTÜR FOTOMETRİSİ')
plt.legend(loc='upper left')
plt.grid(color
         ='white', alpha=0.4, linestyle=':')
plt.tight_layout()
plt.show()