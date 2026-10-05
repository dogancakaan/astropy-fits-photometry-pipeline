# Astropy Tabanlı FITS Verisi Analizi ve Apertür Fotometrisi Pipeline'ı

Bu Python projesi, astrofiziksel veri analizi süreçlerini simüle etmek ve otomatikleştirmek amacıyla **Astropy** kütüphanesini kullanarak sentetik bir FITS (Flexible Image Transport System) dosyası oluşturur, WCS (World Coordinate System) tanımlar ve gelişmiş bir fotometri pipeline'ı çalıştırır.

## Proje Akışı ve Özellikler
1. **Sentetik FITS Üretimi:** Arka plan gürültüsü (`Gaussian` dağılımlı) içeren ve üzerine nokta yayılma fonksiyonu (`PSF - Point Spread Function`) ile yıldızlar eklenmiş sentetik bir gökyüzü görüntüsü üretir.
2. **WCS Entegrasyonu:** Görüntüye Dünya Koordinat Sistemi (`RA`, `Dec`) bilgileri gömer.
3. **Arka Plan Tahmini:** `sigma_clipped_stats` metodunu kullanarak istatistiksel arka plan gürültüsünü ve medyan seviyesini hesaplar.
4. **Apertür Fotometrisi:** Belirlenen hedef yıldız koordinatlarında maskeleme yaparak net akı (flux) hesabı yapar.
5. **Gökyüzü Dönüşümü:** Piksel koordinatlarını `pixel_to_world_values` fonksiyonu ile gökyüzü koordinatlarına (`RA`, `Dec`) çevirir.
6. **Kadir ve Sönümleme Hesabı:** Sıfır noktası (`Zero Point`) kullanarak gözlenen kadiri hesaplar ve yıldızlararası toz sönümlemesi ($A_V = R_V \cdot E(B-V)$) düzeltmesi uygulayarak gerçek (intrinsic) kadiri elde eder.
7. **WCS Destekli Görselleştirme:** `matplotlib` ve `Astropy WCS` projeksiyonunu kullanarak profesyonel bir gökyüzü haritası çıktısı sunar.

## Gereksinimler
Projeyi çalıştırmak için gerekli kütüphaneler:
```bash
pip install numpy matplotlib astropy scipy