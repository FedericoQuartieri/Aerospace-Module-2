# Dati presi da hyStrath (hy2Foam)

hy2Foam e' il codice con cui sono state fatte le figure del paper
(Casseau et al., Aerospace 2016, 3, 34). Il suo repository e' hyStrath:

- URL: https://github.com/vincentcasseau/hyStrath
- commit: 984e3000a5f8 (ramo master, settembre 2026)
- licenza: GNU GPL v3 (file LICENSE.txt del repository)

## File

- `hTCReactionsQK`: copia del file
  `run/hyStrath/hy2Foam/genericCase/constant/chemDicts/hTCReactionsQK`.
  Sono le 19 reazioni dell'aria a 5 specie (15 dissociazioni e 4 scambi,
  irreversibili) con le costanti QK di Scanlon et al., AIAA Journal 53 (2015),
  tabella 1, che il paper cita per la fig. 9 ma non riporta (la tabella 2 del
  paper ha solo la reazione N2 + N2).
  Unita' del file: A in m3/kmol/s, Ta in K.

- `air5_QK.xml` e `air5_QK_back.xml`: le stesse 19 reazioni nel formato dei
  meccanismi di Mutation++, con A convertito in cm3/mol/s (x 1000). Mutation++
  non accetta in un meccanismo due reazioni irreversibili una inversa
  dell'altra (le considera la stessa reazione), quindi i due scambi inversi
  (O2 + N -> NO + O e NO + N -> N2 + O) stanno nel secondo file e `Test-air5`
  somma le velocita' dei due meccanismi. I file sono collegati da
  `../mutation-data/mechanisms/`, cosi' Mutation++ li trova per nome.

Nel resto del progetto, da hyStrath viene solo il valore della sezione d'urto
V-V di N2-O2 (2.667e-19 m2, file `thermo2TModel`), citato nel commento di
`N2O2/Test-N2O2.C`.
