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
hifiasm -o HG002.asm --h1 read1.fq.gz --h2 read2.fq.gz HG002-HiFi.fq.gz
```

- Evaluation
```shell

```

## 2. Haplotype-aware scaffolding
- Remove Haplotype-specific HiC reads
- Scaffold
