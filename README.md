# Praca magisterska - Analiza bezpieczeństwa kryptosystemów opartych na krzywych eliptycznych

Repozytorium zawiera implementacje wszystkich algorytmów wykorzystanych podczas eksperymentów oraz komplet danych wejściowych i otrzymanych wyników.

Implementacje zostały wykonane w języku Python z wykorzystaniem pakietu SageMath. Do generowania wykresów wykorzystano bibliotekę Matplotlib. SageMath udostępnia funkcje służące do wykonywania działań na punktach krzywych eliptycznych, obliczania liczby punktów krzywej oraz wielu innych operacji algebraicznych, dzięki czemu umożliwia wygodną realizację eksperymentów związanych z ECDLP.

## Wymagania

Do uruchomienia dołączonych implementacji wymagane jest zainstalowanie pakietu SageMath oraz biblioteki Matplotlib.

Programy można uruchamiać bezpośrednio poleceniem:

```bash
sage nazwa_pliku.py
```

z poziomu katalogu zawierającego odpowiednie pliki eksperymentu.

## Struktura archiwum

Archiwum zostało podzielone na trzy katalogi odpowiadające przeprowadzonym eksperymentom:

- `bsgs_rho` -- eksperymenty dla algorytmów baby-step giant-step oraz rho Pollarda,
- `mov` -- eksperymenty związane z redukcją MOV,
- `ph` -- eksperymenty dla algorytmu Pohliga--Hellmana.

Każdy katalog zawiera kod źródłowy, dane wejściowe oraz wyniki uzyskane podczas eksperymentów.

## Zawartość katalogów

Każdy z katalogów `bsgs_rho`, `mov` oraz `ph` zawiera następujący zestaw plików:

- `*_collect_data.py` -- generowanie danych testowych oraz przeprowadzanie eksperymentów,
- `*_dataset.csv` -- parametry krzywych eliptycznych i podgrup,
- `*_from_dataset.py` -- przeprowadzanie eksperymentów dla wcześniej przygotowanego zbioru danych,
- `*_results.csv` -- wyniki pomiarów oraz dodatkowe wielkości wyznaczane w trakcie eksperymentów,
- `*_make_plots.py` -- generowanie wykresów i statystyk,
- `*_summary.txt` -- zbiorcze podsumowanie eksperymentu.

## Uwagi dotyczące eksperymentów

Pliki `*_collect_data.py` generują dane testowe, przeprowadzają eksperymenty, a następnie zapisują parametry wygenerowanych krzywych i podgrup do plików `*_dataset.csv` oraz wyniki pomiarów do plików `*_results.csv`.

Pliki `*_from_dataset.py` umożliwiają ponowne przeprowadzenie eksperymentów dla wcześniej przygotowanego zbioru danych zapisanego w plikach `*_dataset.csv`. Otrzymane wyniki zapisywane są do plików `*_results.csv`.

Podczas wykonywania eksperymentów losowane są punkty należące do rozważanych podgrup oraz odpowiadające im logarytmy dyskretne. W związku z tym kolejne uruchomienia programów mogą prowadzić do nieznacznie różnych czasów wykonania algorytmów, nawet przy wykorzystaniu tego samego zbioru krzywych i podgrup.

Pliki `*_make_plots.py` wykorzystują dane zapisane w plikach `*_results.csv` do wygenerowania wykresów oraz statystyk podsumowujących przeprowadzone eksperymenty.