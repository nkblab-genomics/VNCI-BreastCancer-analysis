#!/bin/bash
sample=</sample id>
bam=</bam path>
telseq -H -m -r 150 "$bam" > "$sample".telseq.out
