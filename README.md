# amiR_target_finder

amiR_target_finder identifies minimum sets of artificial micoRNA/s that can target all sequences in a FASTA format reference seqeunce file.

<code>amiR_target_finder [-h] [-min_win MIN_WIN] [-max_win MAX_WIN] [-V]
                          input_file

amiR_target_finder

  Created by Stephen Fletcher on 2016-03-02.
  Copyright 2016 Stephen Fletcher. All rights reserved.

  Licensed under the MIT License 

  Distributed on an "AS IS" basis without warranties
  or conditions of any kind, either express or implied.

USAGE

positional arguments:
  input_file            FASTA file containing reference sequences

optional arguments:
  -h, --help            show this help message and exit
  -min_win MIN_WIN, --min_win MIN_WIN
                        Min. window size to iterate from
  -max_win MAX_WIN, --max_win MAX_WIN
                        Max. window size to iterate to
  -V, --version         show program's version number and exit
scmb-bcarr05m:amiR_target_finder steve$ </code>
