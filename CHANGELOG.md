# Changelog
All notable changes to this project will be documented in this file. See [conventional commits](https://www.conventionalcommits.org/) for commit guidelines.

- - -
## [v0.10.1](https://github.com/yoctoyotta1024/CLEO/compare/9d74ef1f02a0c325e895399a3ea02d8851cc69a1..v0.10.1) - 2024-06-06
#### Bug Fixes
- update yaxt version in CI - ([0e31ed3](https://github.com/yoctoyotta1024/CLEO/commit/0e31ed3ca9e7ef8561b1c266c9ab8681c09527b2)) - clara.bayley
- update yac version in CI - ([0b1a4bd](https://github.com/yoctoyotta1024/CLEO/commit/0b1a4bdd7773066ac4827bdf8acce37417d7f333)) - clara.bayley
#### Refactoring
- update flags and compilers used to install yac - ([24302fc](https://github.com/yoctoyotta1024/CLEO/commit/24302fcc97bd9bf7ea20e02eeca43811eb8609d4)) - clara.bayley
- update yaxt version - ([9d74ef1](https://github.com/yoctoyotta1024/CLEO/commit/9d74ef1f02a0c325e895399a3ea02d8851cc69a1)) - clara.bayley

- - -

## [v0.10.0](https://github.com/yoctoyotta1024/CLEO/compare/7423cef79d566206c38bf00eec41ae2c42f4aeb0..v0.10.0) - 2024-06-06
#### Bug Fixes
- fix LAPACK FLAGS to off for cvode sundials options - ([a1c9a31](https://github.com/yoctoyotta1024/CLEO/commit/a1c9a3144b59504bf1d06e1a3bdfc113e292c503)) - clara.bayley
- more messages for YAC install directory and fix to LAPACK dependency with SUNDIALS - ([09fce52](https://github.com/yoctoyotta1024/CLEO/commit/09fce5277d3e4a856036d69eda4a1f99fce1eee2)) - clara.bayley
- add packages to load if building with YAC - ([87da022](https://github.com/yoctoyotta1024/CLEO/commit/87da022e011fbd773737b3c674c42b12f28138ec)) - clara.bayley
- missing set of enableyac in examples bash scripts - ([bb25192](https://github.com/yoctoyotta1024/CLEO/commit/bb25192c9d0565036e435051421383eef19ebf53)) - clara.bayley
- typo in cmake flag for setting C compiler - ([7c4298a](https://github.com/yoctoyotta1024/CLEO/commit/7c4298a1ba73134edcf9d636f6d3f31c253da125)) - clara.bayley
#### Build system
- added the ENABLE_YAC_COUPLING flag to guard the coupldyn_yac build - ([00792a6](https://github.com/yoctoyotta1024/CLEO/commit/00792a67e66408a15f016c9e42f369cffe6bcdb2)) - wiltonloch
#### Continuous Integration
- updated yac version to v3.2.0_a_p2 - ([b512d28](https://github.com/yoctoyotta1024/CLEO/commit/b512d287466b23690732bc7b3a325b6e3ba3acff)) - wiltonloch
- added divfree2D_yac example as a build workflow step - ([7150405](https://github.com/yoctoyotta1024/CLEO/commit/715040531cc4d532aaca469ef658e745c7c1c04d)) - wiltonloch
#### Documentation
- corrected instruction for requirements to run CLEO with YAC - ([30bef82](https://github.com/yoctoyotta1024/CLEO/commit/30bef824a5ac50841c81743599b7b891db8ed49b)) - clara.bayley
- note on configuring examples to run with yac or not - ([0abb9bf](https://github.com/yoctoyotta1024/CLEO/commit/0abb9bf04100820cbc5ec3e7fbd193a7af6b7f3c)) - clara.bayley
- updated info in packages required for YAC, plus modify changelog to remove annoying message - ([282062c](https://github.com/yoctoyotta1024/CLEO/commit/282062c53b1197dc19c976432f755d8ab02541a5)) - clara.bayley
- how to install YAC info - ([f3703ff](https://github.com/yoctoyotta1024/CLEO/commit/f3703ff16e61870ee1ec5583009e9180fa97aaa8)) - clara.bayley
- how to install YAC info - ([ffffda3](https://github.com/yoctoyotta1024/CLEO/commit/ffffda31972a006798cfcc46d77dfb22378e6350)) - clara.bayley
#### Features
- new files created for easily running yac_3d example - ([94456b1](https://github.com/yoctoyotta1024/CLEO/commit/94456b19c21fbc559a9d5b3d875d377d80a27f46)) - clara.bayley
- option to build YAC added to build CLEO helper scripts - ([2b4a79b](https://github.com/yoctoyotta1024/CLEO/commit/2b4a79b0349a468e454effd5944e2db9dc175c74)) - clara.bayley
- option to build YAC added to build CLEO helper scripts - ([de5e7ec](https://github.com/yoctoyotta1024/CLEO/commit/de5e7ec45880c47a0d279863ae01af0633a6b56d)) - clara.bayley
- new bash script for installing YAC and YAXT on levante - ([af94452](https://github.com/yoctoyotta1024/CLEO/commit/af94452980b76808911154aeb0b0d3f9fa509d65)) - clara.bayley
- new bash script for installing YAC and YAXT on levante - ([dbb9f4d](https://github.com/yoctoyotta1024/CLEO/commit/dbb9f4dfca7797967bfa2b3c1cb68e457f877228)) - clara.bayley
#### Miscellaneous Chores
- renaming - ([54c42d9](https://github.com/yoctoyotta1024/CLEO/commit/54c42d91c14c68decbf0000454de438c0a6a46d6)) - clara.bayley
- rename files - ([7423cef](https://github.com/yoctoyotta1024/CLEO/commit/7423cef79d566206c38bf00eec41ae2c42f4aeb0)) - clara.bayley
#### Refactoring
- add required to find packages for YAC build - ([b552086](https://github.com/yoctoyotta1024/CLEO/commit/b5520860b9422044bf4c337a9986d059044890af)) - clara.bayley
- delete work-in-progress scripts - ([30ff4f5](https://github.com/yoctoyotta1024/CLEO/commit/30ff4f5443b4f84130c987c63e3773225d27465f)) - clara.bayley
- update yac version in installation and minor rename to bash - ([73e2aae](https://github.com/yoctoyotta1024/CLEO/commit/73e2aaea7b4ef17a600fc1405b48a5386a970403)) - clara.bayley
- copy inputfiles python script from fromfile example - ([bac06ab](https://github.com/yoctoyotta1024/CLEO/commit/bac06abb841ac2413bbe0918ced9db5622e6674e)) - clara.bayley
- amend comment - ([77acf3b](https://github.com/yoctoyotta1024/CLEO/commit/77acf3b715eee6f8906982cfd60d2afb61dff975)) - clara.bayley
- update bash scripts to pass in enableyac flag to CLEO build script - ([8f23fb2](https://github.com/yoctoyotta1024/CLEO/commit/8f23fb29ac2e3eb505397867c3d921fc608175c1)) - clara.bayley
- edits to building with YAC scripts to pass enable/disable flag to cmake with/without yac and yaxt root directory - ([b53936b](https://github.com/yoctoyotta1024/CLEO/commit/b53936bef609b5015e7b8a1c57ec5277acdf6942)) - clara.bayley

- - -

## [v0.9.0](https://github.com/yoctoyotta1024/CLEO/compare/46c7b2af9c3ad3d4b1b5ecd451c574f05f5bdbe7..v0.9.0) - 2024-05-26
#### Bug Fixes
- way to include sdmmonitor in microphysics concept - ([4e3c21d](https://github.com/yoctoyotta1024/CLEO/commit/4e3c21d1db98af04b9c4ae427ad9c3ae3b5b2c9f)) - clara.bayley
- way to include sdmmonitor in microphysics concept - ([47891b8](https://github.com/yoctoyotta1024/CLEO/commit/47891b8291b348dea02c97c75bd5ad44464b2830)) - clara.bayley
- add SDMMo template parameter to microphysics concept - ([eb4b868](https://github.com/yoctoyotta1024/CLEO/commit/eb4b868de471f2de7f82d7bfd4a91b0954babb34)) - clara.bayley
#### Build system
- added more general MPI support to use YAC on CXX - ([46c7b2a](https://github.com/yoctoyotta1024/CLEO/commit/46c7b2af9c3ad3d4b1b5ecd451c574f05f5bdbe7)) - wiltonloch
#### Continuous Integration
- added the yac_3d example as a build workflow step - ([1da2dea](https://github.com/yoctoyotta1024/CLEO/commit/1da2dea94f62ab8e50e4cc2471f1d2c19d7f94d0)) - wiltonloch
#### Documentation
- class not struct fix typo - ([80a0e5d](https://github.com/yoctoyotta1024/CLEO/commit/80a0e5de64d6a3e5871eb438526cd4b04c62b693)) - clara.bayley
- better title to page - ([12cdb86](https://github.com/yoctoyotta1024/CLEO/commit/12cdb861b572018bf8ca73dd8e77a02b1aabf964)) - clara.bayley
#### Features
- monitor for condensed mass outputs mass condensed per gbx for each output timestep - ([7b0e644](https://github.com/yoctoyotta1024/CLEO/commit/7b0e64403f8aef31c646008894490ca09b6c6cd4)) - clara.bayley
- use C++ traits and templated functions to determine dtype string for zarr metadata and required xarray metadata given data type - ([ae017fa](https://github.com/yoctoyotta1024/CLEO/commit/ae017fa66dfad677fe96b2c18d3988765c4c166e)) - clara.bayley
- option to use monitor observer in main - ([4cf90a8](https://github.com/yoctoyotta1024/CLEO/commit/4cf90a834a238915c497410291e60acccb9adf29)) - clara.bayley
- placeholder struct for monitoring condensation - ([41a1fda](https://github.com/yoctoyotta1024/CLEO/commit/41a1fda47784034bf401eec0319e93f8ab5af7b3)) - clara.bayley
- struct for combining sdmmonitors - ([69172f3](https://github.com/yoctoyotta1024/CLEO/commit/69172f31a93a4607a9937afab3c2bbe8e0105b21)) - clara.bayley
- new struct for detecting variables within sdm subtimestepped processes - ([152fa89](https://github.com/yoctoyotta1024/CLEO/commit/152fa895d3c21c77183a4a229af265f223329537)) - clara.bayley
- new function constraint on observers concept - ([62a1a3c](https://github.com/yoctoyotta1024/CLEO/commit/62a1a3c3543418a6e0d744ce3666e57904095012)) - clara.bayley
#### Miscellaneous Chores
- correct typo in comment - ([85f8072](https://github.com/yoctoyotta1024/CLEO/commit/85f8072041452b12087686dd1d13b7e4939bf79e)) - clara.bayley
- correct comment and formatting - ([b46a0fc](https://github.com/yoctoyotta1024/CLEO/commit/b46a0fc137d3a975a4af21c925d3eec6a5aa09f8)) - clara.bayley
- delete excess notes - ([03cfd09](https://github.com/yoctoyotta1024/CLEO/commit/03cfd099097f8e8aff153868eb101ca0f8af58b9)) - clara.bayley
- rename const tstep observer - ([b40e55b](https://github.com/yoctoyotta1024/CLEO/commit/b40e55b7f0e3d65ddb0956495936234c81666d28)) - clara.bayley
- rename - ([4c289ee](https://github.com/yoctoyotta1024/CLEO/commit/4c289eefab5b2b2d994399c615b0867abecdae52)) - clara.bayley
- rename - ([a0a21d1](https://github.com/yoctoyotta1024/CLEO/commit/a0a21d14978987a3dfbaf4fd71587e5bd6e09af3)) - clara.bayley
- rename file - ([d159b90](https://github.com/yoctoyotta1024/CLEO/commit/d159b90888d56a80fbce53a52215d8eb429e9d41)) - clara.bayley
#### Refactoring
- work to test monitor microphysics in condensation - ([b65bacd](https://github.com/yoctoyotta1024/CLEO/commit/b65bacdce397e7a9436d3c4b9fbe0f4b51a6879e)) - clara.bayley
- revert to commit e3041238475c3d4ac848ce8bb34b6f9f95ddbf07 and then move sdmmonitor to superdrops library - ([9336266](https://github.com/yoctoyotta1024/CLEO/commit/933626698b94d875b3bdd8bfeade888e7d3770f4)) - clara.bayley
- Revert "refactor: attempt to include monitor in condensation" - ([147cc9d](https://github.com/yoctoyotta1024/CLEO/commit/147cc9d5bc1a78a28d31876737c447933db39266)) - clara.bayley
- attempt to include monitor in condensation - ([7ac5712](https://github.com/yoctoyotta1024/CLEO/commit/7ac5712e3a70fc8184fd1e54b73a110b62eb69a6)) - clara.bayley
- sdm microphysics concept requires sdm monitor - ([02af1d2](https://github.com/yoctoyotta1024/CLEO/commit/02af1d2581fa7a3710f5b9343e121aafcd43d4e7)) - clara.bayley
- use explicit template to aid cuda compilation - ([8dde10e](https://github.com/yoctoyotta1024/CLEO/commit/8dde10ee7485f81e78993b258a61d42e4059e5ae)) - clara.bayley
- make sdm monitor observer generic - ([5d0b1a0](https://github.com/yoctoyotta1024/CLEO/commit/5d0b1a01b1ed28d88acbffd80eda084705354a32)) - clara.bayley
- use buffer view type in monitor - ([bf61b09](https://github.com/yoctoyotta1024/CLEO/commit/bf61b091722c554a3fe6ec8266393286723c2a46)) - clara.bayley
- sketch for outline of condensation observer - ([824a501](https://github.com/yoctoyotta1024/CLEO/commit/824a50133bec3568008ac9778ec0683c27b42594)) - clara.bayley
- move sdmmonitors into subdirectory - ([b570d21](https://github.com/yoctoyotta1024/CLEO/commit/b570d217827d370a0337b44bbef7863d57192de1)) - clara.bayley
- convert SDMMonitor into concept - ([e54c75f](https://github.com/yoctoyotta1024/CLEO/commit/e54c75fac2196f0ce5d2788be4ef04c3264797e3)) - clara.bayley
- convert SDMMonitor into concept - ([83decb3](https://github.com/yoctoyotta1024/CLEO/commit/83decb3106f18ee4d649cdb88406bbda90a26011)) - clara.bayley
- convert SDMMonitor into concept - ([32d146d](https://github.com/yoctoyotta1024/CLEO/commit/32d146d65ed6cdc2f6945582a2dda8d5881b69bd)) - clara.bayley
- delete redundant files - ([d2ab7f0](https://github.com/yoctoyotta1024/CLEO/commit/d2ab7f06a9ecf050ed15606b0e9f95b0a11f1674)) - clara.bayley
- add empty func to satisfy new const step observer - ([85cc173](https://github.com/yoctoyotta1024/CLEO/commit/85cc1734caf24d81be5ff9b517f52f6cea1e3c79)) - clara.bayley
- const step observer also calls sdm_step func - ([21dbebc](https://github.com/yoctoyotta1024/CLEO/commit/21dbebcda7089edf0f6894fc81ffe5e913f6e9f3)) - clara.bayley
- seperate files for const tstep observer and rename struct - ([ebfe99e](https://github.com/yoctoyotta1024/CLEO/commit/ebfe99e083f57bbe7b551949203810ce0cfb35e6)) - clara.bayley

- - -

## [v0.8.2](https://github.com/yoctoyotta1024/CLEO/compare/afa7b400d9b4e214239df293784e1695824a5378..v0.8.2) - 2024-05-24
#### Bug Fixes
- corrected build and execution for yac_3d example - ([5d3d147](https://github.com/yoctoyotta1024/CLEO/commit/5d3d14733ed4f37665425c17c57a4be194b33114)) - wiltonloch
#### Build system
- moved mpi and fyaml requirements from coupldyn_yac to findyac - ([60239f3](https://github.com/yoctoyotta1024/CLEO/commit/60239f3273a0dd5e1ab8d3b6d2153915c81274d3)) - wiltonloch
#### Refactoring
- move random number generator into collisions - ([afa7b40](https://github.com/yoctoyotta1024/CLEO/commit/afa7b400d9b4e214239df293784e1695824a5378)) - clara.bayley

- - -

## [v0.8.1](https://github.com/yoctoyotta1024/CLEO/compare/a4c22eb26a906faef8d91f140fc229d7025fd537..v0.8.1) - 2024-05-08
#### Bug Fixes
- docs corrected message abotu ruff settings - ([3e55894](https://github.com/yoctoyotta1024/CLEO/commit/3e558940718e085da1dfa849785dce1e250ed2ff)) - clara.bayley
- value error is index error - ([cf9a8da](https://github.com/yoctoyotta1024/CLEO/commit/cf9a8da63a2fd06904255a0685071309c07cb4b6)) - clara.bayley
- add Wextra to debugging compiler flags - ([3d1eacb](https://github.com/yoctoyotta1024/CLEO/commit/3d1eacbd9f49a4128b26ad041261d6b5959b950a)) - clara.bayley
- wrong type of error except - ([fa353c3](https://github.com/yoctoyotta1024/CLEO/commit/fa353c319b3cb93a024b1f2d7872b5556af069da)) - clara.bayley
- Key not Value Errors for dicts and datasets - ([4f0e2e4](https://github.com/yoctoyotta1024/CLEO/commit/4f0e2e49805d6eb883afa1b14cbbbe0a4735da70)) - clara.bayley
- Key not Value Errors for dicts and datasets - ([274be91](https://github.com/yoctoyotta1024/CLEO/commit/274be913378060604d1e750fd0c3900e3d925ffd)) - clara.bayley
- ignore Module level import not at top of file error in ruff linting - ([440f8d9](https://github.com/yoctoyotta1024/CLEO/commit/440f8d996779a042cdce49e39d652c30df2ece41)) - clara.bayley
#### Documentation
- include info about ruff in docs - ([a4c22eb](https://github.com/yoctoyotta1024/CLEO/commit/a4c22eb26a906faef8d91f140fc229d7025fd537)) - clara.bayley
#### Miscellaneous Chores
- format and linting Python with ruff - ([5b09491](https://github.com/yoctoyotta1024/CLEO/commit/5b09491d04ce0bc5b4820d1aac8732cf9c6e1f0c)) - clara.bayley
- format and linting Python with ruff - ([dc4233d](https://github.com/yoctoyotta1024/CLEO/commit/dc4233db06602f7f7315009b7627785ba0c7a957)) - clara.bayley
- format and linting Python with ruff - ([4ddb02f](https://github.com/yoctoyotta1024/CLEO/commit/4ddb02f8ef70896dcb4d6828ae7ea19d14eb647b)) - clara.bayley

- - -

## [v0.8.0](https://github.com/yoctoyotta1024/CLEO/compare/cb5447f387e109abffbf5b21b4522423bd6e0540..v0.8.0) - 2024-05-07
#### Features
- use ruff in pre-commit for python linting and formatting - ([cb5447f](https://github.com/yoctoyotta1024/CLEO/commit/cb5447f387e109abffbf5b21b4522423bd6e0540)) - clara.bayley

- - -

## [v0.7.0](https://github.com/yoctoyotta1024/CLEO/compare/f765ebb289fe4a4359f6e67c1c070dca0432733a..v0.7.0) - 2024-05-06
#### Bug Fixes
- exclude top value in interval - ([ef288a7](https://github.com/yoctoyotta1024/CLEO/commit/ef288a72b5ad7e64f9048a8175d8f245ca99ad55)) - clara.bayley
#### Features
- add supers bc is gpu compatible but has strange warning... - ([796e8a0](https://github.com/yoctoyotta1024/CLEO/commit/796e8a063769b8ca154829d3a417da1318741785)) - clara.bayley
#### Miscellaneous Chores
- testing with longer run works on cpus and gpus - ([ceb3220](https://github.com/yoctoyotta1024/CLEO/commit/ceb322076c314cbaca1fc2ecfaaba500627f4696)) - clara.bayley
- correct includes - ([597c4b0](https://github.com/yoctoyotta1024/CLEO/commit/597c4b0bc14244dd49adf684bd048f5473ae5f7c)) - clara.bayley
#### Refactoring
- setup eurec4a1D example with different boundary conditions - ([48ed274](https://github.com/yoctoyotta1024/CLEO/commit/48ed27471e79c7ff7db7ffb71dc049809ce98392)) - clara.bayley
- better comments explaining encapsulation - ([fe761db](https://github.com/yoctoyotta1024/CLEO/commit/fe761dbe1718218a3957a6c388e56ec5363bbdef)) - clara.bayley
- encapsulate motion functions into struct within movement - ([74ccef7](https://github.com/yoctoyotta1024/CLEO/commit/74ccef7dee104f3a3dad055bbe22d3467e08d48d)) - clara.bayley
- move constructor to cpp - ([0bc5863](https://github.com/yoctoyotta1024/CLEO/commit/0bc5863885612c14616b2d20a4556090d11d67d3)) - clara.bayley
- use std rand gen not kokkos - ([a15494f](https://github.com/yoctoyotta1024/CLEO/commit/a15494f39f58e0be0f4649ddc39947e7c67c2fde)) - clara.bayley
- remove supers returns gbxindexes view for adding supers - ([ab21b49](https://github.com/yoctoyotta1024/CLEO/commit/ab21b49be8dedace691bac25a75f85779c36b9ed)) - clara.bayley
- parallel verison of removing supers - ([d7c3afa](https://github.com/yoctoyotta1024/CLEO/commit/d7c3afa577eb80769b869a925b3b72d6030cc1a7)) - clara.bayley
- move funcs out of class - ([3d3fa6e](https://github.com/yoctoyotta1024/CLEO/commit/3d3fa6e48c5cf45b59196e76dd014786319f0fbf)) - clara.bayley
- use add supers bcs - ([f765ebb](https://github.com/yoctoyotta1024/CLEO/commit/f765ebb289fe4a4359f6e67c1c070dca0432733a)) - clara.bayley

- - -

## [v0.6.1](https://github.com/yoctoyotta1024/CLEO/compare/84e58e0b90123b3f15cf31233d5a9509ea2cfb64..v0.6.1) - 2024-05-03
#### Bug Fixes
- typo falsely renaming function reverted - ([418b387](https://github.com/yoctoyotta1024/CLEO/commit/418b387c9c728ee9d1f7e3aad2bf08579bf85509)) - clara.bayley
- move gbxvol and gbxarea maps to host from device - ([3d5959a](https://github.com/yoctoyotta1024/CLEO/commit/3d5959ab0305b8838d93c6e5f10e9d766808768c)) - clara.bayley
- move gbxvol and gbxarea maps to host from device - ([84e58e0](https://github.com/yoctoyotta1024/CLEO/commit/84e58e0b90123b3f15cf31233d5a9509ea2cfb64)) - clara.bayley

- - -

## [v0.6.0](https://github.com/yoctoyotta1024/CLEO/compare/e2dfb97594856e94a6aad8d9d47fb322835c46dd..v0.6.0) - 2024-05-03
#### Bug Fixes
- update examples with renamed config var - ([3346d33](https://github.com/yoctoyotta1024/CLEO/commit/3346d339d2995e2120216572bbf595b64f3cf585)) - clara.bayley
- update examples with renamed config var - ([e01fb44](https://github.com/yoctoyotta1024/CLEO/commit/e01fb4481db828ba43222e6200518b781f022ffd)) - clara.bayley
#### Features
- seperate script for making yac1 input files - ([f49cd97](https://github.com/yoctoyotta1024/CLEO/commit/f49cd9742dd9febf0cd30aef23ee4086fc972c74)) - clara.bayley
- seperate script for making yac1 input files - ([e2dfb97](https://github.com/yoctoyotta1024/CLEO/commit/e2dfb97594856e94a6aad8d9d47fb322835c46dd)) - clara.bayley
#### Refactoring
- split initial conditions generation from divfree2d run and plotting - ([b02c3aa](https://github.com/yoctoyotta1024/CLEO/commit/b02c3aaeb1b82c53543bdeed67f78bf0f9478b61)) - clara.bayley
- seperate input file generation from running fromfile model and plotting of model output - ([d25323b](https://github.com/yoctoyotta1024/CLEO/commit/d25323b5bc21d070edc19f290df258c7ecab1e23)) - clara.bayley

- - -

## [v0.5.1](https://github.com/yoctoyotta1024/CLEO/compare/82b968815ee53a229fd869315bcf3076b56cc30e..v0.5.1) - 2024-05-03
#### Bug Fixes
- missing macros for GPU compilation added - ([82b9688](https://github.com/yoctoyotta1024/CLEO/commit/82b968815ee53a229fd869315bcf3076b56cc30e)) - clara.bayley

- - -

## [v0.5.0](https://github.com/yoctoyotta1024/CLEO/compare/e59ea81a340ac66c6e9e121dab03b02a22d76ff1..v0.5.0) - 2024-05-02
#### Bug Fixes
- restore changelog since v0.1.0 - ([7281aa1](https://github.com/yoctoyotta1024/CLEO/commit/7281aa15085968fd22d4915a344453e53d2bd29c)) - clara.bayley
#### Features
- corrected cocogitto workflow wtith git push to add changelog changes after bump - ([e59ea81](https://github.com/yoctoyotta1024/CLEO/commit/e59ea81a340ac66c6e9e121dab03b02a22d76ff1)) - clara.bayley
#### Refactoring
- formating of commands in post bump hooks - ([59a0357](https://github.com/yoctoyotta1024/CLEO/commit/59a035753cfe91395b98a4807b4b61d061cb2068)) - clara.bayley

- - -

## [v0.4.0](https://github.com/yoctoyotta1024/CLEO/compare/v0.3.2..v0.4.0) - 2024-05-02
#### Features
- specify repo for release - ([8f215d2](https://github.com/yoctoyotta1024/CLEO/commit/8f215d2ee7fac864cd3435bcc52bbafe65b7a59b)) - clara.bayley

- - -

## [v0.3.2](https://github.com/yoctoyotta1024/CLEO/compare/v0.3.1..v0.3.2) - 2024-05-02
#### Bug Fixes
- bad token - ([1d8e3a5](https://github.com/yoctoyotta1024/CLEO/commit/1d8e3a5d9ddac61719a792576e16efc785ef3131)) - clara.bayley
- minor changes to CI to try to make changelog release work - ([c48e2c0](https://github.com/yoctoyotta1024/CLEO/commit/c48e2c0169a36ba26e377a0c158155dcfd88ee1f)) - clara.bayley

- - -

## [v0.3.1](https://github.com/yoctoyotta1024/CLEO/compare/v0.3.0..v0.3.1) - 2024-05-02
#### Bug Fixes
- test cocogitto working with changelog updates on version bump - ([bfa66d8](https://github.com/yoctoyotta1024/CLEO/commit/bfa66d85becf4725416edf73826d9c56f911ec8d)) - clara.bayley
#### Refactoring
- don't specify github token - ([9f69343](https://github.com/yoctoyotta1024/CLEO/commit/9f69343c6e0bb5056d567384b909c86b4431cc48)) - clara.bayley
- delete explicity repo and add explicit token - ([a4753b8](https://github.com/yoctoyotta1024/CLEO/commit/a4753b84492ea396126bd25b5be0d8667fea9c54)) - clara.bayley
- change repo name - ([a4aafb5](https://github.com/yoctoyotta1024/CLEO/commit/a4aafb56dd700940c1c4cb099e448f73cee2c9d8)) - clara.bayley
- attempt to fix changelog publication with repository spec - ([7695d91](https://github.com/yoctoyotta1024/CLEO/commit/7695d91ae3d9b89aca4cee6c354e2ee05c3b481f)) - clara.bayley
- correct changelog since v0.1.0 - ([5e0f721](https://github.com/yoctoyotta1024/CLEO/commit/5e0f721544d6a5c960dbc340b3b9af4b2e93e573)) - clara.bayley
- correct changelog since v0.1.0 - ([0aae8e9](https://github.com/yoctoyotta1024/CLEO/commit/0aae8e93723098b8d6e03637d2bfebb0e74c350f)) - clara.bayley

- - -

## [v0.3.0](https://github.com/yoctoyotta1024/CLEO/compare/v0.2.1..v0.3.0) - 2024-05-02
#### Bug Fixes
- use remote template instead of full_hash to match cog.toml - ([8529c75](https://github.com/yoctoyotta1024/CLEO/commit/8529c750462fd50fdce2874b9e2fc601fb70f073)) - clara.bayley
#### Features
- pre and post bump messages in cocogitto CI - ([71594cc](https://github.com/yoctoyotta1024/CLEO/commit/71594cccea48b7748a0351d1a860ca6d71b31ec1)) - clara.bayley

- - -

## [v0.2.1](https://github.com/yoctoyotta1024/CLEO/compare/v0.2.0..v0.2.1) - 2024-05-02
#### Bug Fixes
- changelog from v0.1.0 to v0.2.0 manually updated - ([d9052b8](https://github.com/yoctoyotta1024/CLEO/commit/d9052b8e9a08fef017accd357f567c3ef00c3dc0)) - clara.bayley

- - -

## [v0.2.0](https://github.com/yoctoyotta1024/CLEO/compare/c358c9069e701b42f1d33f5cdf6ad672e6aa03a5..v0.2.0) - 2024-05-02
#### Features
- changelog file added and release of new version included to github CI - ([6fdc113](https://github.com/yoctoyotta1024/CLEO/commit/6fdc1130cc5e1b79f602e5fe4f30e67f2fff5fe5)) - clara.bayley
- cocogitto configuration file - ([2990d63](https://github.com/yoctoyotta1024/CLEO/commit/2990d6321dd84895dfa57cc06653eb3996a53a35)) - clara.bayley
- new file for cocogitto in CI on pushes - ([c358c90](https://github.com/yoctoyotta1024/CLEO/commit/c358c9069e701b42f1d33f5cdf6ad672e6aa03a5)) - clara.bayley

- - -

This changelog was generated by [cocogitto](https://github.com/oknozor/cocogitto).
