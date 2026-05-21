# Simulasi Mekanisme Transpor Membran Sel

Tugas Besar IF3211 Komputasi Domain-Spesifik
Institut Teknologi Bandung

---

## Deskripsi Program

Program ini mensimulasikan dan membandingkan tiga mekanisme transpor molekul melintasi membran sel:

1. **Difusi Sederhana** menggunakan Hukum Fick untuk memodelkan pergerakan pasif molekul searah gradien konsentrasi.
2. **Difusi Terfasilitasi** menggunakan kinetika Michaelis-Menten untuk memodelkan transpor yang dimediasi protein, termasuk efek saturasi pada konsentrasi tinggi.
3. **Transpor Aktif** menggunakan sistem ODE untuk mensimulasikan pompa Na+/K+ ATPase yang mempertahankan gradien ion melawan arah gradien konsentrasi dengan menggunakan energi ATP.

Selain simulasi utama, program juga melakukan analisis sensitivitas untuk mengidentifikasi parameter mana yang paling berpengaruh terhadap perilaku sistem, serta perbandingan kondisi normal versus abnormal (sel sehat vs sel dengan ATP rendah akibat iskemia).

---

## Tujuan

- Mengimplementasikan model komputasi dari proses biologis nyata menggunakan Python.
- Membandingkan perilaku ketiga mekanisme transpor secara kuantitatif dan visual.
- Mengidentifikasi parameter kritis melalui analisis sensitivitas lokal (koefisien elastisitas).
- Menunjukkan dampak kondisi patologis (kekurangan ATP) terhadap homeostasis ion sel.

---

## Fitur

- Simulasi difusi sederhana dengan solusi analitik eksak (stabil secara numerik).
- Simulasi difusi terfasilitasi dengan solver ODE scipy (RK45).
- Simulasi pompa Na+/K+ dengan parameter biologis yang konsisten secara fisiologis.
- Analisis sensitivitas untuk semua parameter di ketiga model.
- Perbandingan kondisi normal vs abnormal (ATP normal vs ATP rendah).
- Empat grafik output berkualitas tinggi yang tersimpan otomatis sebagai file PNG.
- Ringkasan hasil simulasi dan ranking parameter tercetak ke konsol.

---

## Struktur File

```
membrane_transport/
|-- simple_diffusion.py       Simulasi difusi sederhana (Hukum Fick)
|-- facilitated_diffusion.py  Simulasi difusi terfasilitasi (Michaelis-Menten)
|-- active_transport.py       Simulasi pompa Na+/K+ ATPase (sistem ODE)
|-- sensitivity_analysis.py   Analisis sensitivitas untuk ketiga model
|-- visualization.py          Pembuatan dan penyimpanan semua grafik
|-- main.py                   Entry point utama program
|-- outputs/
    |-- plot1_concentration_vs_time.png
    |-- plot2_flux_vs_concentration.png
    |-- plot3_sensitivity_analysis.png
    |-- plot4_normal_vs_abnormal.png
```

---

## Persyaratan

- Python 3.10 atau lebih baru
- numpy
- scipy
- matplotlib

Install dependensi dengan perintah:

```
pip install numpy scipy matplotlib
```

---

## Cara Penggunaan

Jalankan program utama dari direktori `membrane_transport/`:

```
cd membrane_transport
python main.py
```

Program akan:
1. Menjalankan ketiga simulasi transpor membran.
2. Melakukan analisis sensitivitas untuk semua model.
3. Mencetak ringkasan hasil dan ranking parameter ke konsol.
4. Menyimpan empat grafik PNG ke folder `outputs/`.

Untuk menjalankan atau mengeksplorasi modul secara individual:

```
python simple_diffusion.py
python facilitated_diffusion.py
python active_transport.py
python sensitivity_analysis.py
python visualization.py
```

---

## Output Grafik

| File | Deskripsi |
|------|-----------|
| plot1_concentration_vs_time.png | Konsentrasi vs waktu untuk ketiga mekanisme |
| plot2_flux_vs_concentration.png | Perbandingan fluks vs konsentrasi: difusi sederhana vs terfasilitasi (menunjukkan saturasi) |
| plot3_sensitivity_analysis.png | Grafik batang ranking parameter berdasarkan indeks sensitivitas |
| plot4_normal_vs_abnormal.png | Overlay kondisi normal vs ATP rendah pada pompa Na+/K+ |

---

## Parameter Biologis Baseline

| Parameter | Nilai | Keterangan |
|-----------|-------|------------|
| D (koef. difusi) | 1e-9 m2/s | Molekul polar kecil di bilayer lipid |
| L (tebal membran) | 7e-9 m | Bilayer fosfolipid tipikal |
| C_awal difusi sederhana | 100 mM vs 10 mM | Gradien konsentrasi antar kompartemen |
| Jmax (difusi terfasilitasi) | 1e-5 mol/m2/s | Orde besaran transporter GLUT1 |
| Km (difusi terfasilitasi) | 5 mM | Konstanta afinitas transporter |
| Na+ intrasel (normal) | 15 mM | Nilai fisiologis mamalia |
| Na+ ekstrasel | 145 mM | Nilai fisiologis mamalia |
| K+ intrasel (normal) | 140 mM | Nilai fisiologis mamalia |
| K+ ekstrasel | 5 mM | Nilai fisiologis mamalia |
| ATP (normal) | 5 mM | Sel sehat |
| ATP (abnormal) | 0.5 mM | Kondisi iskemia/kekurangan energi |

Sumber: Alberts et al., Molecular Biology of the Cell, 6th ed.; Koeppen & Stanton, Berne & Levy Physiology, 6th ed.

---

## Pembagian Tugas

| No | Nama | NIM | Tugas Utama | Modul | Keterangan |
|----|------|-----|-------------|-------|------------|
| 1 | [Nama 1] | [NIM 1] | Implementasi difusi sederhana dan difusi terfasilitasi | simple_diffusion.py, facilitated_diffusion.py | |
| 2 | [Nama 2] | [NIM 2] | Implementasi transpor aktif (pompa Na+/K+) | active_transport.py | |
| 3 | [Nama 3] | [NIM 3] | Analisis sensitivitas dan validasi parameter biologis | sensitivity_analysis.py | |
| 4 | [Nama 4] | [NIM 4] | Visualisasi, integrasi modul, dan dokumentasi | visualization.py, main.py, README.md | |

---

## Referensi

- Alberts, B. et al. (2014). Molecular Biology of the Cell, 6th ed. Garland Science.
- Koeppen, B.M. & Stanton, B.A. (2009). Berne & Levy Physiology, 6th ed. Mosby Elsevier.
- Lodish, H. et al. (2016). Molecular Cell Biology, 8th ed. W.H. Freeman.
- Lauger, P. (1991). Electrogenic Ion Pumps. Sinauer Associates.
- Skou, J.C. (1957). The influence of some cations on an adenosine triphosphatase from peripheral nerves. Biochim. Biophys. Acta, 23, 394-401.
- Saltelli, A. et al. (2008). Global Sensitivity Analysis: The Primer. Wiley.
