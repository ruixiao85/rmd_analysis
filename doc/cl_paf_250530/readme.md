# ChengLei AF 2024 data


```ps1

pandoc -f docx -t markdown --extract-media=. ".\paper\Figures.docx" -o figures.md

pandoc -f docx -t markdown --extract-media=. ".\paper\Cover letter.docx" -o "Cover letter.md"
pandoc -f docx -t markdown --extract-media=. ".\paper\Precision HF prevention in PAF.docx" -o "Precision HF prevention in PAF.md"
pandoc -f docx -t markdown --extract-media=. ".\paper\Supplementary tables.docx" -o "Supplementary tables.md"
pandoc -f docx -t markdown --extract-media=. ".\paper\Tables and figures.docx" -o "Tables and figures.md"


```