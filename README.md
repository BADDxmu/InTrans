## TransIntegrator ： a pipeline for integrative transcript library construction
## Integrative transcript library for <em>Branchiostoma floridae</em> (www.bio-add.org/InTrans/)
### Copyright (C) 2017 ZhiLiang Ji (appo@xmu.edu.cn)
## RequirementT
This software is suitable for all unix-like system with python(version 2.7.7) installed.<br>
`One python module was required before usage : configparser3.5.0.`<br><br>
Moreover, three already published softwares should be correctly installed in advance, and
make sure they had been add to your system environment variables. The three softwares are:<br>
(1) IDBA (version 1.1.1) <br>
(2) CD-HIT (version 4.5.4) <br>
(3) CAP3 (version 12/21/07) <br>
of course, for softwares mentioned above, other version is allowed. However, the pipeline operated
stably with the recommended version. <br>
## Installation Guide
Simply installed by extracting the software package
## Usage
In the package folder you extracted, there are three files and one derectory : `"InTrans.py", "run.cfg", "__init__.py" and "test_data"`<br>
(1) "InTrans.py" is the software `executed file`<br>
(2) "run.cfg" is the `configure file`, which contains a series of important parameters. For correctly running with your data, you set the right parameter value in "run.cfg" file. `Detail of these parameters is writed in "run.cfg"`, or if you confused, please see the corresponding software manual.
### Warnning
(1) `the default maximun read length of IDBA  is 128 bp`, if your read is longer than that, you should change
the vaue of 128 to longer one (e.g. 250) in "xx/idba-xxx/src/sequence/short_sequence.h" : <br>
"static const uint32_t kMaxShortSequence = 128;"<br>
-><br>
"static const uint32_t kMaxShortSequence = 250;" <br>
(2) correspondingly, you should also change the default kmer unit to bigger one(e.g. 8) in "xx/idba-xxxsrc/basic/kmer.h":<br>
"static const uint32_t kNumUint64 = 4;"<br>
-><br>
"static const uint32_t kNumUint64 = 8;"<br>
(3) recompile IDBA after modification to make new read length and kmer working <br>
## Running
If individual parameter value had been set in "run.cfg" file, then run the pipeline with: <br>
$ `python InTrans.py run.cfg`<br>
For example, you can make a test running with datas in `"test_data"`:<br>
(1) run without heterogeneous data, corresponding configure file is `run_fq.cfg`:<br>
$ `cd ./test_data/`<br>
$ `python ../InTrans.py run_fq.cfg` <br>
(2) run with heterogeneous data, corresponding configure file is `run_fq_heterogen.cfg`:<br>
$ `cd ./test_data/`<br>
$ `python ../InTrans.py run_fq_heterogen.cfg` <br>
## Output
Two folders and one log file were generated after the program runs out: <br>
(1) "output" folder <br>
    contains the final transcript file, which in fasta format.<br>
(2) "temp_output" folder<br>
    contains the temporary file during running, include output of IDBA, CD-HIT, and CAP3. <br>
