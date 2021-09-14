# phased-assembly-check

## 1. Haplotype-resolved ctg

- Trio
```shell
yak count -b37 -t16 -o pat.yak <(cat pat_1.fq.gz pat_2.fq.gz) <(cat pat_1.fq.gz pat_2.fq.gz)
yak count -b37 -t16 -o mat.yak <(cat mat_1.fq.gz mat_2.fq.gz) <(cat mat_1.fq.gz mat_2.fq.gz)
hifiasm -o HG002.asm -t32 -1 pat.yak -2 mat.yak HG002-HiFi.fa.gz
```
- HiC
```shell
hifiasm -o HG002.asm -t32 --h1 read1.fq.gz --h2 read2.fq.gz HG002-HiFi.fq.gz
```

- Evaluation
```shell
## yak (https://github.com/lh3/yak)

yak trioeval pat.yak mat.yak HG002.asm.trio.hap1.p_ctg.fa > HG002.asm.hap1.trioeval
yak trioeval pat.yak mat.yak HG002.asm.trio.hap2.p_ctg.fa > HG002.asm.hap2.trioeval

## Merqury.FK (https://github.com/thegenemyers/MERQURY.FK)

```

## 2. Haplotype-aware scaffolding
- Remove Haplotype-specific HiC reads (https://www.nature.com/articles/s41586-021-03451-0)

```shell
meryl-lookup -memory 2 -exclude -mers pat.meryl -sequence $read1 -sequence2 $read2 -r2
mat.R2.fastq.gz | pigz -c > mat.R1.fastq.gz
meryl-lookup -memory 2 -exclude -mers mat.meryl -sequence $read1 -sequence2 $read2 -r2
pat.R2.fastq.gz | pigz -c > pat.R1.fastq.gz 

Using https://github.com/marbl/merqury/blob/master/trio/exclude_reads.sh
```
- Scaffold two haplotype together and check
```
1. juicer + 3d-dna / AllHiC
2. Loading into JuiceBox for visualization

- False duplication

(Examples are based on the hifiasm 0.15.4 for high het plant genome with trio data)

    - Fig1
    > ![Fig1](images/Fig1.png)
    > |     seqName        |     #matKmer    |     #patKmer    |     #pat-pat    |     #pat-mat    |     #mat-pat    |     #mat-mat    |     seqLen     |
    > |--------------------|-----------------|-----------------|-----------------|-----------------|-----------------|-----------------|----------------|
    > |     h1tg000041l    |          160    |        10644    |           77    |           83    |           82    |        10561    |     6152553    |
    > |     h1tg000053l    |        50105    |           62    |        50059    |           45    |           45    |           17    |     5159680    |
    > |     h1tg000050l    |        17945    |           27    |        17924    |           20    |           20    |            7    |     2082689    | 
    Conclusion: Remove the h1tg000041l

    - Fig2
> 
>|     seqName        |     #patKmer    |     #patKmer    |     #pat-pat    |     #pat-mat    |     #mat-pat    |     #mat-mat    |     seqLen     |
>|--------------------|-----------------|-----------------|-----------------|-----------------|-----------------|-----------------|----------------|
>|     h1tg000015l    |           87    |          164    |           51    |           35    |           36    |          128    |     2770227    |
>|     h1tg000078l    |           42    |           66    |           24    |           18    |           18    |           47    |     1476449    |
>|     h1tg000028l    |         7532    |           87    |         7478    |           53    |           53    |           34    |     3799565    |

    - Fig3
    - Fig4
    - Fig5

- Misplaced haplotype
If you find the `hap1` ctg have more strong interaction with the hap2, and `trioeval` support it, you need manually move `this hap1` ctg to hap2



```



