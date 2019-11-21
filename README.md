Инструкции для запуска ДЗ 1 и ДЗ 2:

    **ДЗ 1:** 
    python -m alignment.aligner -seq1 GGTTGACTA -seq2 TGTTACGGTA -alignment -gap -2 — **глобальное выравнивание простое**
    python -m alignment.test_aligner — **тест**
    
    **ДЗ 2:**
    python -m alignment.aligner -seq1 GGTTGACTA -seq2 TGTTACGGTA -alignment -gap -2 -local — **локальное выравнивание**
    python -m alignment.aligner_with_affine_gaps -seq1 GGTTGACTA -seq2 TGTTACGGTA -alignment -gap_start -2 -gap_continued -0.5 — **глобальное с аффинными штрафами**
    
    python -m alignment.random_nucleotides — **генератор**
