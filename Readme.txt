// Catalin Constantin Usurelu
// 341 C1

Aplicatia se poate rula deschizand fisierul index.html intr-un browser cu suport HTML 5 (de preferabil o versiune recenta de Chrome).
De asemenea aplicatia se poate accesa la adresa http://128.199.38.110:8080/ (in caz ca serverul este inactiv folositi varianta locala).

Interfata aplicatiei este intuitiva. Butonul choose file permite alegera fisierului cu training samples (ca cele furnizate de dumneavoastra).
Campurile de completat sunt parametrii pentru functiile de aflat minimumEmbeddingDimension si de prediction (numele sunt corespunzatoare cu
cele din curs). Daca nu sunt completate vor fi folosite valorile default inscrise in ele. Observatie: graficul initial (cel afisat
imediat ce intram in aplicatie este initializat cu o functie care nu are nici-o legatura cu tema, este doar de initializare).

Pe ultimul rand este scrisa mean square error + numele fisierului pentru care am obtinut rezultatul (in caz ca am uitat pentru ce calculam).

Prin plasarea cursorului pe grafic se afiseaza valoarea pentru predictie/original in acel punct (se poate da si zoom prin scroll).

Mentionez ca aplicatia merge cam greu pentru date de intrare mari ( > 2000 de sample-uri), deci s-ar putea sa tina 5-10 secunde o predictie.
In orice caz, chiar daca sta mai mult, in final va furniza rezultatul. De asemenea, tot pentru date mari, si graficul va afisa mai greu (cursor/scroll).

De asemenea, am atasat un fisier cu date de intrarea care mi-a fost foarte util pentru debugging si sunt convins ca le-ar fi si altor daca l-ar detine (este
o functie care reprezinta o oscilatie in rezonanta). Stiu ca nu este ceva care contine noise si nici nu e un semnal nestationar, dar scopul este
de a exemplifica modul de predictie si de a ajuta la debugging (este mai vizual rezultatul).

!Observatie - in caz ca apar erori sau graficul este incomplet vor trebui introduse alte date de intrare (marire angle si micsorare/marire 
nearest neighbours, deoare pentru anumite combinatii este posibil sa nu avem suficienti nearest neighbours si sa nu putem face predictia).


Fisierul de interes este: js/Nearest Neighbour Predictor.js