Rast provodnika – Python alat

Ovaj projekat predstavlja interaktivni Python alat za proračun rastojanja između dva provodnika u elektroenergetskim vodovima.
Aplikacija koristi Tkinter grafički interfejs za unos parametara i prikazuje 2D i 3D grafikone sa rezultatima proračuna.

Funkcionalnosti:

-Unos parametara za dva provodnika putem Tkinter GUI-ja

-Proračun oblika lančanice, ugiba i ugla otklona

-Automatsko generisanje 2D i 3D grafikona

-Izračunavanje potrebnih i stvarnih rastojanja između provodnika

-Identifikacija kritičnih tačaka i provera usklađenosti

-Automatski izvoz rezultata u Excel (.xlsx)

Korišćene biblioteke:

-numpy

-sympy

-scipy

-matplotlib

-tkinter

-pandas

-colorama

Instalacija potrebnih biblioteka:

pip install numpy sympy scipy matplotlib pandas colorama


(tkinter je deo standardne Python instalacije.)

Pokretanje aplikacije:

Klonirajte ili preuzmite repozitorijum:

git clone https://github.com/korisnik/rast-provodnika.git
cd rast-provodnika


Pokrenite program:

python "rast provodnika_@ ver 0.py"


Unesite potrebne podatke (stubovi, dužine, parametar lančanice, pritisak vetra, tip provodnika itd.)

Nakon potvrde unosa, program automatski:

-prikazuje 2D i 3D grafikone provodnika,

-izračunava potrebna i stvarna rastojanja,

-identifikuje kritične tačke,

-čuva rezultate u Excel fajl Rastojanje provodnika VP.xlsx.

Namena:

Ovaj alat je razvijen kao pomoć elektroenergetskim inženjerima za mehaničke proračune nadzemnih vodova i verifikaciju rastojanja između provodnika.
