#!/bin/bash

mtx="/local/jrev_suite-symmetrics/*.mtx.sorted"

mail.pl  kkourt@cslab.ece.ntua.gr "Here we go on $(hostname)" ...

scripts/crs32_run_clone.sh          $mtx | tee ~/J.rev/csr32-$(hostname)-$(iddate)
mail.pl  kkourt@cslab.ece.ntua.gr "CSR done on $(hostname)" ...

scripts/csrdu_run_clone.sh          $mtx | tee ~/J.rev/csrdu-$(hostname)-$(iddate)
mail.pl  kkourt@cslab.ece.ntua.gr "CSRDU done on $(hostname)" ...

scripts/csrvi_all_run_clone.sh      $mtx | tee ~/J.rev/csrvi-$(hostname)-$(iddate)
mail.pl  kkourt@cslab.ece.ntua.gr "CSRVI done on $(hostname)" ...

scripts/csrdu_vi_all_run_clone.sh   $mtx | tee ~/J.rev/csrduvi-$(hostname)-$(iddate)
mail.pl  kkourt@cslab.ece.ntua.gr "CSRDUVI done on $(hostname)" ...

scripts/bcsr32_run_clone.sh         $mtx | tee ~/J.rev/bcsr32-$(hostname)-$(iddate)
mail.pl  kkourt@cslab.ece.ntua.gr "BSCR done on $(hostname)" ...
