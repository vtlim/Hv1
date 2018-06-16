/* vim: set filetype=markdown : */



# taut 1

* VANILLA
  * Loading
    ```
    vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/structural/analyzeDCD.tcl -args ../../taut1_18629-19_s211.psf 1 ../../taut1_18629-19_s211.pdb ../../npt01-fwd.dcd
    ```
  * Analysis
    ```
    calc_rmsf_hv1 rmsf_hv1
    count_wat_z waters-in-zrange.dat protein,and,resid,159,and,name,CA protein,and,resid,118,and,name,CA
    ```

* FEP
  * Loading
    ```
    vmdt -e /beegfs/DATA/mobley/limvt/hv1/04_fep/analysis/structural/analyzeDCD.tcl -args ../../00_main/18629-19_R211S.psf 1 ../../00_main/18629-19_R211S.pdb ../../FEP_F/lambda_01/alchemy01.dcd ../../FEP_F/lambda_02/alchemy02.dcd ../../FEP_F/lambda_03/alchemy03.dcd ../../FEP_F/lambda_04/alchemy04.dcd ../../FEP_F/lambda_05/alchemy05.dcd ../../FEP_F/lambda_06/alchemy06.dcd ../../FEP_F/lambda_07/alchemy07.dcd ../../FEP_F/lambda_08/alchemy08.dcd ../../FEP_F/lambda_09/alchemy09.dcd ../../FEP_F/lambda_10/alchemy10.dcd ../../FEP_F/lambda_11/alchemy11.dcd ../../FEP_F/lambda_12/alchemy12.dcd ../../FEP_F/lambda_13/alchemy13.dcd ../../FEP_F/lambda_14/alchemy14.dcd ../../FEP_F/lambda_15/alchemy15.dcd ../../FEP_F/lambda_16/alchemy16.dcd ../../FEP_F/lambda_17/alchemy17.dcd ../../FEP_F/lambda_18/alchemy18.dcd ../../FEP_F/lambda_19/alchemy19.dcd ../../FEP_F/lambda_20/alchemy20.dcd
    ```
