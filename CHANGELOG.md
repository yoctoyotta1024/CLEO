# Changelog
All notable changes to this project will be documented in this file. See [conventional commits](https://www.conventionalcommits.org/) for commit guidelines.

- - -
## [v0.44.3](https://github.com/yoctoyotta1024/CLEO/compare/88166164197cf06f0cda37ea463dfe49cf0385e1..v0.44.3) - 2025-06-27
#### Bug Fixes
- no longer need finalize in example - ([77a5029](https://github.com/yoctoyotta1024/CLEO/commit/77a5029fa44b1751316d0a7fe49a9d932f2076af)) - clara.bayley
- construct arrays by reference not copy - ([8816616](https://github.com/yoctoyotta1024/CLEO/commit/88166164197cf06f0cda37ea463dfe49cf0385e1)) - clara.bayley

- - -

## [v0.44.2](https://github.com/yoctoyotta1024/CLEO/compare/c1ab081fb2a53556b85a0389de8af438c2250161..v0.44.2) - 2025-06-26
#### Bug Fixes
- prevent multiple kokkos init/finalize in python bindings with at exit call - ([34ac7af](https://github.com/yoctoyotta1024/CLEO/commit/34ac7af6f26a25d11b35c09b8a47b58545cfdff3)) - clara.bayley
#### Documentation
- update YAC requirements information - ([a64cda0](https://github.com/yoctoyotta1024/CLEO/commit/a64cda0cac4420e1d12acb4bbf4bc7e779eab069)) - clara.bayley
#### Refactoring
- python bindings example uses new generalised SDMMethods - ([e5ac1dc](https://github.com/yoctoyotta1024/CLEO/commit/e5ac1dcc5323861fe145ca7f8fb6245e0c27955e)) - clara.bayley
- binding for combined null and condensation microphysical process - ([363ed08](https://github.com/yoctoyotta1024/CLEO/commit/363ed084f9533d25f390952de82bf912ccdca1ca)) - clara.bayley
- generalise sdm methods bindings - ([bf82ca3](https://github.com/yoctoyotta1024/CLEO/commit/bf82ca344dcf17c66e770e251f53b02a46863c2e)) - clara.bayley
- python bound function to return a combination of null and condensation microphysics - ([f3474dd](https://github.com/yoctoyotta1024/CLEO/commit/f3474dd9a9621eaf505734fb54f029d5bfa53a8a)) - clara.bayley
- add sdm methods bindings for condensation-only - ([8890be0](https://github.com/yoctoyotta1024/CLEO/commit/8890be0fc6b8baf049a2fe35823518cfea7d0552)) - clara.bayley
- special case for maximum interval of microphysics - ([753b32e](https://github.com/yoctoyotta1024/CLEO/commit/753b32eb0373b52ff56ccf57f23f6a39c250b40c)) - clara.bayley
- add guard on kokkos initialise - ([9753575](https://github.com/yoctoyotta1024/CLEO/commit/97535751ffdf14bffa032545d4ec41da25c39851)) - clara.bayley
- module not spack to load openmpi and use newer intel compiler - ([c1ab081](https://github.com/yoctoyotta1024/CLEO/commit/c1ab081fb2a53556b85a0389de8af438c2250161)) - clara.bayley

- - -

## [v0.44.1](https://github.com/yoctoyotta1024/CLEO/compare/636db598409bf012ff1c50b506efc8dd221e7a2e..v0.44.1) - 2025-06-23
#### Bug Fixes
- zarr library does not depend on cartesian decomposition - ([a6eb12a](https://github.com/yoctoyotta1024/CLEO/commit/a6eb12ac1ee0ce4ff32587c32f5ae8b066d1374f)) - clara.bayley
- zarr lib depends on cartesiandomain - ([dc18395](https://github.com/yoctoyotta1024/CLEO/commit/dc18395f9a38f6097242390169b12fc1fcde5126)) - clara.bayley
- use gcc compiler in run examples - ([3499eb9](https://github.com/yoctoyotta1024/CLEO/commit/3499eb9e35c37c86a698f0549807823c1562485a)) - clara.bayley
- bubble requires two mpi tasks in SLURM - ([f702f7f](https://github.com/yoctoyotta1024/CLEO/commit/f702f7f77cae03300495f7de0b3c8476768b4e33)) - clara.bayley
- Fixed a compiler warning related error (fallthrough error) - ([636db59](https://github.com/yoctoyotta1024/CLEO/commit/636db598409bf012ff1c50b506efc8dd221e7a2e)) - clara.bayley
#### Continuous Integration
- all examples not need yaxt/yac in build - ([284af75](https://github.com/yoctoyotta1024/CLEO/commit/284af7520cc30f07343ffcd569c8f67030037ccd)) - clara.bayley
#### Miscellaneous Chores
- rename Dataset -> Simple or Collective Dataset - ([cac0dbf](https://github.com/yoctoyotta1024/CLEO/commit/cac0dbf323c7e38128c2d9cb3a31dc190dbe3ca3)) - clara.bayley
- add todo - ([feb4b49](https://github.com/yoctoyotta1024/CLEO/commit/feb4b49d34f50218748e2dc2a7ba24409df3922a)) - clara.bayley
#### Performance Improvements
- remove unwanted comment - ([bc00384](https://github.com/yoctoyotta1024/CLEO/commit/bc00384b73e8eb22b36c65b9772349ece71a82a3)) - clara.bayley
#### Refactoring
- generalise operator to combine two CollectDataForDataset types - ([e80fde8](https://github.com/yoctoyotta1024/CLEO/commit/e80fde842154fb3212e22f3e3d7ec50ca9ca978b)) - clara.bayley
- make examples compatible with templated dataset - ([ec19973](https://github.com/yoctoyotta1024/CLEO/commit/ec19973dde412b08d5bdbd0ca424d0087b692bdf)) - clara.bayley
- template over dataset - ([4bc7cdd](https://github.com/yoctoyotta1024/CLEO/commit/4bc7cdde00fd59386293130f74f0f569209f3f29)) - clara.bayley
- remove guard on collective dataset - ([69c70f3](https://github.com/yoctoyotta1024/CLEO/commit/69c70f353d646e698cc50c93c496e0c422513d40)) - clara.bayley
- fix sbatch tasks and nthreads for bubble test case - ([ee350a4](https://github.com/yoctoyotta1024/CLEO/commit/ee350a4113865d3d3acc2a8a7e9b5e9527fa8282)) - clara.bayley
- make pycleo compatible with new configuration library - ([a343c3f](https://github.com/yoctoyotta1024/CLEO/commit/a343c3f751c18d20d9d9f98d68d40cf562141c28)) - clara.bayley
- move configuration related files into seperate library from initialisation - ([c663c30](https://github.com/yoctoyotta1024/CLEO/commit/c663c30bd1eba2819148e497fbec45d14b7e8713)) - clara.bayley
- change order of cmake building - ([2dbca64](https://github.com/yoctoyotta1024/CLEO/commit/2dbca64753fc9324a53831bfe20b44bded688f61)) - clara.bayley
- generalise yac installation to allow intel compiler on Levante - ([711714c](https://github.com/yoctoyotta1024/CLEO/commit/711714c8bc72fa4f8365b0b8fdb1b4ed4bd88211)) - clara.bayley
- delete redundant enableyacpython flag - ([72e4577](https://github.com/yoctoyotta1024/CLEO/commit/72e4577034a2d3c8614549a791a97797f05bd86d)) - clara.bayley
- make building yac an essential requirement of CLEO build scripts - ([b61c1bb](https://github.com/yoctoyotta1024/CLEO/commit/b61c1bb19e3060ff083178358d2fbf87b70e437a)) - clara.bayley
- enableyac -> enable_yacpython flag renaming - ([8883c60](https://github.com/yoctoyotta1024/CLEO/commit/8883c6085e06b46b71f14e4d158051c90f98a487)) - clara.bayley
- long time for python bindings example - ([dbbfc71](https://github.com/yoctoyotta1024/CLEO/commit/dbbfc71d9447fa9ba2849c51b9a5fdc7158a4524)) - clara.bayley

- - -

## [v0.44.0](https://github.com/yoctoyotta1024/CLEO/compare/3b38b2c1c6665568cb10eb08d48130a06796d2d6..v0.44.0) - 2025-06-12
#### Bug Fixes
- coupldyn_numpy is submodule of pycleo - ([a0c42bb](https://github.com/yoctoyotta1024/CLEO/commit/a0c42bb4f85597b6342132367d72e55455882bc8)) - clara.bayley
#### Features
- python class for thermodynamics of example - ([9dc56f3](https://github.com/yoctoyotta1024/CLEO/commit/9dc56f3ded58362a2fec8d22b30933332e1ca869)) - clara.bayley
- new library for numpy arrays coupled dynamics - ([439289e](https://github.com/yoctoyotta1024/CLEO/commit/439289e2c1b727d26e2a534296b32a5ed8242384)) - clara.bayley
#### Miscellaneous Chores
- renaming config conflicting variables - ([42d4cff](https://github.com/yoctoyotta1024/CLEO/commit/42d4cff787356cf4b1706c18f2d8a2fe3a19b155)) - clara.bayley
#### Performance Improvements
- no inline in macro - ([cc7ae0d](https://github.com/yoctoyotta1024/CLEO/commit/cc7ae0d0d8d9a4f6d25403a33a66704c001a487b)) - clara.bayley
#### Refactoring
- use coupldyn_numpy submodule in python_bindings example - ([de975a9](https://github.com/yoctoyotta1024/CLEO/commit/de975a95d5224b9bf3c7a2d081cd0fa6aad0f48c)) - clara.bayley
- use thermodynamics in python bindings example - ([91e12e5](https://github.com/yoctoyotta1024/CLEO/commit/91e12e599f94057e3314d5c0096f90840e9c3557)) - clara.bayley
- parallelise numpy comms send/receive - ([bfeff8e](https://github.com/yoctoyotta1024/CLEO/commit/bfeff8ead50997d83862a805ba26c88b3bf7f1ca)) - clara.bayley
- bindings for time to model timestep conversion - ([01ca4f2](https://github.com/yoctoyotta1024/CLEO/commit/01ca4f28ba70ce6aefe93b6daf14c75b955e97ff)) - clara.bayley
- reduce number of gridboxes - ([038f858](https://github.com/yoctoyotta1024/CLEO/commit/038f858a1e502341873a3cb674aa050f9084377c)) - clara.bayley
- delete unused thermofiles coupled dynamics - ([0751f80](https://github.com/yoctoyotta1024/CLEO/commit/0751f80a546d217b4f27028604939ecdcfacf940)) - clara.bayley
- use cleo config struct in intialisation of cleo not python config - ([b7fdde4](https://github.com/yoctoyotta1024/CLEO/commit/b7fdde49ee44e2bddbb9dd93808111cceae8dd5b)) - clara.bayley
- include more getters in config bindings - ([3b38b2c](https://github.com/yoctoyotta1024/CLEO/commit/3b38b2c1c6665568cb10eb08d48130a06796d2d6)) - clara.bayley

- - -

## [v0.43.0](https://github.com/yoctoyotta1024/CLEO/compare/0a032e280dcee5934cfa2430884569131fa27875..v0.43.0) - 2025-06-12
#### Features
- first bindings for gridboxes - ([6480b59](https://github.com/yoctoyotta1024/CLEO/commit/6480b591fe6dac629b2df42cf4dd00ef1d30a27e)) - clara.bayley
- first bindings for superdroplets - ([085d86c](https://github.com/yoctoyotta1024/CLEO/commit/085d86c956a9d82142b3698c769ae252961a886a)) - clara.bayley
- new python bindings for initialisation/configuration - ([90bd3fc](https://github.com/yoctoyotta1024/CLEO/commit/90bd3fcf17463d2ada5180f13fb2de92a6f24188)) - clara.bayley
- first bindings for boundary conditions and movesupersindomain - ([21b595a](https://github.com/yoctoyotta1024/CLEO/commit/21b595ad739d398ab20a761e4f84344df7b183e5)) - clara.bayley
- first bindings for cartesian transport - ([59ce95c](https://github.com/yoctoyotta1024/CLEO/commit/59ce95cb12df1a017d30dd86f99788b1a1c569e6)) - clara.bayley
- create first bindings for motion and microphysics - ([9eb782a](https://github.com/yoctoyotta1024/CLEO/commit/9eb782a496c48c8e10d06f1b8c4c578708a2bb36)) - clara.bayley
- create first binding for cartesian maps - ([be74b04](https://github.com/yoctoyotta1024/CLEO/commit/be74b04b42c3613c1121412d6a4041adfc16b9c6)) - clara.bayley
- create first bindings for observers - ([57cdc4e](https://github.com/yoctoyotta1024/CLEO/commit/57cdc4ebb248d028ce41c63aae08869640e40fdc)) - clara.bayley
- create first bindings for SDMMethods - ([0a032e2](https://github.com/yoctoyotta1024/CLEO/commit/0a032e280dcee5934cfa2430884569131fa27875)) - clara.bayley
#### Miscellaneous Chores
- renaming and formatting - ([b9a81c5](https://github.com/yoctoyotta1024/CLEO/commit/b9a81c5bb6d8520a617fe6b885050ccdd9aa9956)) - clara.bayley
- add placeholder notes on rest of objects to create - ([0dfff04](https://github.com/yoctoyotta1024/CLEO/commit/0dfff0498f1a9ceaf39cd9d11dfd1953174674b4)) - clara.bayley
#### Refactoring
- create SDMMethods in python bindings example - ([85c2a42](https://github.com/yoctoyotta1024/CLEO/commit/85c2a42ecf6628cfe5eff64dcd014bf98a05826a)) - clara.bayley
- add bindings to timestep functions called by timestep_cleo - ([223ab55](https://github.com/yoctoyotta1024/CLEO/commit/223ab554f3d865f15f99006ca84fb128b33e8707)) - clara.bayley
- add sdm timestepping routines to sdmmethods bindings - ([2ed1589](https://github.com/yoctoyotta1024/CLEO/commit/2ed158993c5a7c102d51f8b90e6263fdd3f21bde)) - clara.bayley
- add access to gbxmaps from SDMMethods - ([c3a4fd8](https://github.com/yoctoyotta1024/CLEO/commit/c3a4fd8db086e5cf173c98822af5114f0731d04a)) - clara.bayley
- add function signature to cartesian maps - ([85375e5](https://github.com/yoctoyotta1024/CLEO/commit/85375e51d527ca2a325d9add554653c594d033c7)) - clara.bayley
- kokkos init takes config - ([8792764](https://github.com/yoctoyotta1024/CLEO/commit/87927643ac889caac8577053755ae33053f5a9c4)) - clara.bayley
- kokkos init and finalise via pycleo - ([03279f6](https://github.com/yoctoyotta1024/CLEO/commit/03279f657ffc90dfa3a1002f78104749dbba36b5)) - clara.bayley

- - -

## [v0.42.0](https://github.com/yoctoyotta1024/CLEO/compare/17ffe2acf119fde7ac55ed2f87132a0429c34a37..v0.42.0) - 2025-06-11
#### Bug Fixes
- add archive library destination to libs targets - ([0b07ff4](https://github.com/yoctoyotta1024/CLEO/commit/0b07ff4bbcb52007242f60d7b70905a3fc913b3a)) - clara.bayley
#### Documentation
- add pybind11 dependency description - ([34514d0](https://github.com/yoctoyotta1024/CLEO/commit/34514d0c154b8d9533447fec2f56af3eabf5acb0)) - clara.bayley
#### Features
- add new example for testing python bindings - ([bfce784](https://github.com/yoctoyotta1024/CLEO/commit/bfce78478f1819f0f95d822f4623b74dd62c2ebc)) - clara.bayley
- new flag to not build python bindings - ([b8f3ea7](https://github.com/yoctoyotta1024/CLEO/commit/b8f3ea7f237339ad1388e4aa9dffc1bd43f25498)) - clara.bayley
- new library for CLEO's python bindings - ([732c77c](https://github.com/yoctoyotta1024/CLEO/commit/732c77cd83e543bc3e28ef77176a6966a0c83511)) - clara.bayley
#### Performance Improvements
- split long lines over multiple - ([2142890](https://github.com/yoctoyotta1024/CLEO/commit/21428902ea59b72a83cef8d9e33f4cb601be2fe9)) - clara.bayley
#### Refactoring
- add option to specify python version for python bindings via cmake - ([4223ea0](https://github.com/yoctoyotta1024/CLEO/commit/4223ea075c8ad9cbd40e5820fd8276985b5547e5)) - clara.bayley
- move next_couplstep function into sdmmethods - ([17ffe2a](https://github.com/yoctoyotta1024/CLEO/commit/17ffe2acf119fde7ac55ed2f87132a0429c34a37)) - clara.bayley

- - -

## [v0.41.2](https://github.com/yoctoyotta1024/CLEO/compare/33b8848188cdacccfb958aac4a40b00e594cc5f2..v0.41.2) - 2025-06-10
#### Bug Fixes
- fix order of gridbox sizes set in cartesian decompositon - ([33b8848](https://github.com/yoctoyotta1024/CLEO/commit/33b8848188cdacccfb958aac4a40b00e594cc5f2)) - clara.bayley

- - -

## [v0.41.1](https://github.com/yoctoyotta1024/CLEO/compare/5afbfe2fba4c4d031642e7608d40021bf61a274d..v0.41.1) - 2025-06-04
#### Bug Fixes
- plotting of 2d motion - ([3cb3537](https://github.com/yoctoyotta1024/CLEO/commit/3cb353765acda16faeb60ae300546b6ef2470c5b)) - clara.bayley
- plotting of ragged data when data not available - ([b89f674](https://github.com/yoctoyotta1024/CLEO/commit/b89f67456cc6bf980c146718342c5fce31c0773e)) - clara.bayley
#### Refactoring
- add option to detach time from superdroplets - ([e6d2d67](https://github.com/yoctoyotta1024/CLEO/commit/e6d2d67194319ea8a50c18afe1fc64603308065b)) - clara.bayley
- make examples compatible with new api - ([5afbfe2](https://github.com/yoctoyotta1024/CLEO/commit/5afbfe2fba4c4d031642e7608d40021bf61a274d)) - clara.bayley

- - -

## [v0.41.0](https://github.com/yoctoyotta1024/CLEO/compare/1bc7564218e7da23f6fee626c839d1893aabc1e1..v0.41.0) - 2025-06-04
#### Features
- new superdrops module for handing ragged superdroplet arrays - ([9f0b2c8](https://github.com/yoctoyotta1024/CLEO/commit/9f0b2c8ebe5571e5abd8ba3653244931daae7b1b)) - clara.bayley
#### Miscellaneous Chores
- formatting - ([edb0e4d](https://github.com/yoctoyotta1024/CLEO/commit/edb0e4d92f3a0bda4e1ee6ad061450db6c68e51b)) - clara.bayley
#### Refactoring
- add function for selecting specific times of superdrop data - ([6882963](https://github.com/yoctoyotta1024/CLEO/commit/688296384651d72815566c965c0f8d3b716644a3)) - clara.bayley
- add option to attach time to superdroplets - ([f49d1f5](https://github.com/yoctoyotta1024/CLEO/commit/f49d1f5c39d49e29a5d54ffda04d027d1db736b7)) - clara.bayley
- add superdrop sampling functions - ([d1cabd2](https://github.com/yoctoyotta1024/CLEO/commit/d1cabd23f93aae313f1bce5d239c67b669483d77)) - clara.bayley
- add depreciation warnings - ([289e3fe](https://github.com/yoctoyotta1024/CLEO/commit/289e3feba2263b43b3ae28cc89b7b5f0da9f6e8d)) - clara.bayley
- move ensembzarr out of pySD into examples - ([7efb7c4](https://github.com/yoctoyotta1024/CLEO/commit/7efb7c417a237eb7d688a19dd42de02ce14637f6)) - clara.bayley
- better expression of supersdata class - ([a92d079](https://github.com/yoctoyotta1024/CLEO/commit/a92d0790ca14ea0c94c7a92f374cf88003f44227)) - clara.bayley
- a test setup with all superdrop observers - ([1bc7564](https://github.com/yoctoyotta1024/CLEO/commit/1bc7564218e7da23f6fee626c839d1893aabc1e1)) - clara.bayley

- - -

## [v0.40.0](https://github.com/yoctoyotta1024/CLEO/compare/29a1591c5c7142307475c5587b028c240e5c86dc..v0.40.0) - 2025-05-29
#### Bug Fixes
- fix superdroplet plotting when domain is smaller than 100km - ([0d11ffe](https://github.com/yoctoyotta1024/CLEO/commit/0d11ffe90d643e7a3d4f38168f2c4e61e86242a2)) - clara.bayley
#### Features
- new script for plotting bubble - ([b2d57d1](https://github.com/yoctoyotta1024/CLEO/commit/b2d57d10c5c05d16cf4c07f0ef0e8c0eb4bb0e25)) - clara.bayley
#### Refactoring
- export yacyaxt root if running example with yac - ([e0c1959](https://github.com/yoctoyotta1024/CLEO/commit/e0c1959cc10d3005e4daa868c6d03c88486cf967)) - clara.bayley
- use gcc compiler when enabled yac in example - ([b096fa3](https://github.com/yoctoyotta1024/CLEO/commit/b096fa3051ddfc362f825553a72146a8c6d85582)) - clara.bayley
- remove unnecesary parts of tmp_working_bubble script - ([8595efd](https://github.com/yoctoyotta1024/CLEO/commit/8595efd3e15f35d6c9bedbdf26c81aa0f8d791f9)) - clara.bayley
- use cleoenv python to run bubble - ([6c645ed](https://github.com/yoctoyotta1024/CLEO/commit/6c645ed7dcdef54ef6062f89820a74574ec59c92)) - clara.bayley
- update bubble script to run and plot bubble - ([a3d31d0](https://github.com/yoctoyotta1024/CLEO/commit/a3d31d07de86071b601eeec1d294db9d9e7fc173)) - clara.bayley
- read config from yac_bubble_data_reader - ([05b24fa](https://github.com/yoctoyotta1024/CLEO/commit/05b24fa52709549aae90f3e981128b1770ac342b)) - clara.bayley
- move icon yac init parameters into yaml config file - ([66080ef](https://github.com/yoctoyotta1024/CLEO/commit/66080ef6d616f1cb29fe6c282127cd2f2602d631)) - clara.bayley
- add alternative bubble plot - ([d33e7ff](https://github.com/yoctoyotta1024/CLEO/commit/d33e7ffb2de0087db3e034919a49b06eedc3d2a9)) - clara.bayley
- add more superdroplet attribute observers - ([4c66d8e](https://github.com/yoctoyotta1024/CLEO/commit/4c66d8e79adbcce3e2bbf4e0dbcaee53933f01ce)) - clara.bayley
- adapt cleo domain for bubble to crop inside central portion of icon domain - ([770d369](https://github.com/yoctoyotta1024/CLEO/commit/770d369dfd602188db66515917d97ef84ce99031)) - clara.bayley
- adapt cleo domain for bubble to match entire icon domain - ([5056765](https://github.com/yoctoyotta1024/CLEO/commit/5056765048a40d75327e41ad3de5398af4bc966b)) - clara.bayley
- add yacyaxt root export to compile/run cleo bash in case yac enabled - ([29a1591](https://github.com/yoctoyotta1024/CLEO/commit/29a1591c5c7142307475c5587b028c240e5c86dc)) - clara.bayley

- - -

## [v0.39.7](https://github.com/yoctoyotta1024/CLEO/compare/dbec7a4c7ea667b5241b1016a3c29f9fb9727bcd..v0.39.7) - 2025-04-19
#### Bug Fixes
- fix major bug in calculation of collision probability - ([dbec7a4](https://github.com/yoctoyotta1024/CLEO/commit/dbec7a4c7ea667b5241b1016a3c29f9fb9727bcd)) - clara.bayley

- - -

## [v0.39.6](https://github.com/yoctoyotta1024/CLEO/compare/ea79f2a57af61c360e3af46091a5971b70c0b668..v0.39.6) - 2025-04-17
#### Bug Fixes
- fix major bug in calculation of collision probability - ([ea79f2a](https://github.com/yoctoyotta1024/CLEO/commit/ea79f2a57af61c360e3af46091a5971b70c0b668)) - clara.bayley
#### Documentation
- remove no longer relevant statement - ([eb2bce7](https://github.com/yoctoyotta1024/CLEO/commit/eb2bce7a06d014f3705bb44336ada7016c440f03)) - clara.bayley
#### Refactoring
- initial conditions as in Shima et al. 2009 - ([ef16652](https://github.com/yoctoyotta1024/CLEO/commit/ef166525726af75bfaa89bcc170d3ec0c1b8d579)) - clara.bayley

- - -

## [v0.39.5](https://github.com/yoctoyotta1024/CLEO/compare/736fbcbe88438e642b358a38915a28c69970a683..v0.39.5) - 2025-04-17
#### Bug Fixes
- error in Long 1974 collision efficieny git status - ([3643b33](https://github.com/yoctoyotta1024/CLEO/commit/3643b333c3e5f298562b06dc66f69e600024056d)) - clara.bayley
#### Refactoring
- return lines not axes from figure - ([736fbcb](https://github.com/yoctoyotta1024/CLEO/commit/736fbcbe88438e642b358a38915a28c69970a683)) - clara.bayley

- - -

## [v0.39.4](https://github.com/yoctoyotta1024/CLEO/compare/25acc999a08780bbc2ee02a395bbfa156f50bb18..v0.39.4) - 2025-04-16
#### Bug Fixes
- 525da8f8c8b9c0732b97175b1b4351b5bb7ec276 need mpi at runtime when using yac - ([4d11a8a](https://github.com/yoctoyotta1024/CLEO/commit/4d11a8ac1b4f0a33780920e4ae718d39ad2d67c8)) - clara.bayley
#### Refactoring
- revert 37c296a4fe937281f7dd95526ed76d19edcbadb7 and add requirements - ([f7b513b](https://github.com/yoctoyotta1024/CLEO/commit/f7b513b2eff96d8ef74f3fa69f90ba2b16b3aa13)) - clara.bayley
- use python given as argument to create YAC python bindings - ([25acc99](https://github.com/yoctoyotta1024/CLEO/commit/25acc999a08780bbc2ee02a395bbfa156f50bb18)) - clara.bayley

- - -

## [v0.39.3](https://github.com/yoctoyotta1024/CLEO/compare/a9dd722d8b1b24b980a851e555368e0fe4ce80be..v0.39.3) - 2025-04-16
#### Bug Fixes
- need mpi4py in using yac at runtime - ([525da8f](https://github.com/yoctoyotta1024/CLEO/commit/525da8f8c8b9c0732b97175b1b4351b5bb7ec276)) - clara.bayley
#### Refactoring
- remove spack packages from yac runtime settings - ([17a5dc5](https://github.com/yoctoyotta1024/CLEO/commit/17a5dc521a10ef5693a0a15b3d0b4fdcd9a1be8d)) - clara.bayley
- lower python version to match than used to make python bindings on levante - ([37c296a](https://github.com/yoctoyotta1024/CLEO/commit/37c296a4fe937281f7dd95526ed76d19edcbadb7)) - clara.bayley
- use levante_packages bash in yac installation - ([a9dd722](https://github.com/yoctoyotta1024/CLEO/commit/a9dd722d8b1b24b980a851e555368e0fe4ce80be)) - clara.bayley

- - -

## [v0.39.2](https://github.com/yoctoyotta1024/CLEO/compare/ff2accbdbc00d2884d9597ba8702896319b69945..v0.39.2) - 2025-04-11
#### Bug Fixes
- update acyronym - ([ff2accb](https://github.com/yoctoyotta1024/CLEO/commit/ff2accbdbc00d2884d9597ba8702896319b69945)) - clara.bayley

- - -

## [v0.39.1](https://github.com/yoctoyotta1024/CLEO/compare/f5c01a097ad62b1a5544bc19b1931dc95b0af24a..v0.39.1) - 2025-04-11
#### Bug Fixes
- use relative path from current script in bash directories - ([c85b4e9](https://github.com/yoctoyotta1024/CLEO/commit/c85b4e98ec11ed07478e215c438c7d6ec8b30626)) - clara.bayley
#### Refactoring
- allow no executables to be compiled - ([1167c0e](https://github.com/yoctoyotta1024/CLEO/commit/1167c0eb9cedc3502bd2acd631480d8bf732a667)) - clara.bayley
- remove mail user from SLURM jobs - ([f5c01a0](https://github.com/yoctoyotta1024/CLEO/commit/f5c01a097ad62b1a5544bc19b1931dc95b0af24a)) - clara.bayley

- - -

## [v0.39.0](https://github.com/yoctoyotta1024/CLEO/compare/57b2f11b3c1a44bd30936477c463a28061aec4bf..v0.39.0) - 2025-03-26
#### Bug Fixes
- remove conflicting black and ruff python formatting - ([7adb597](https://github.com/yoctoyotta1024/CLEO/commit/7adb5976a5aac9c174e042984acb9668b1ad1970)) - clara.bayley
#### Features
- split thermodynamics generation into thermo and wind fields seperately - ([57b2f11](https://github.com/yoctoyotta1024/CLEO/commit/57b2f11b3c1a44bd30936477c463a28061aec4bf)) - clara.bayley
#### Refactoring
- use threads for constthermo2d and in speedtest builds - ([c0f4043](https://github.com/yoctoyotta1024/CLEO/commit/c0f40435c2ceaca1516aeea35efe6e9b81abe552)) - clara.bayley
- make examples compatible with thermogen changes - ([909f471](https://github.com/yoctoyotta1024/CLEO/commit/909f471e4a18772c2b0f685c5517e7f02e607fce)) - clara.bayley

- - -

## [v0.38.3](https://github.com/yoctoyotta1024/CLEO/compare/077501a990f4a7a06e0755a05bedfda095529784..v0.38.3) - 2025-03-25
#### Bug Fixes
- fail-safe way to use modules in bash scripts - ([077501a](https://github.com/yoctoyotta1024/CLEO/commit/077501a990f4a7a06e0755a05bedfda095529784)) - clara.bayley

- - -

## [v0.38.2](https://github.com/yoctoyotta1024/CLEO/compare/b17c906d030eabac61421e979bb8d4c105603444..v0.38.2) - 2025-03-24
#### Bug Fixes
- typo - ([b17c906](https://github.com/yoctoyotta1024/CLEO/commit/b17c906d030eabac61421e979bb8d4c105603444)) - clara.bayley
#### Documentation
- update info on examples - ([4588d2c](https://github.com/yoctoyotta1024/CLEO/commit/4588d2ce88d25d7c28791c1f3c7e7b1350f94651)) - clara.bayley
- include more information in the quickstart - ([bbf0cfd](https://github.com/yoctoyotta1024/CLEO/commit/bbf0cfdab9126b425719ce0fc75248a624765ea0)) - clara.bayley
- update requirements - ([cac5301](https://github.com/yoctoyotta1024/CLEO/commit/cac5301147c81ad783640d7f6508164471f4ecfa)) - clara.bayley
- add title - ([7ceb73a](https://github.com/yoctoyotta1024/CLEO/commit/7ceb73ac63a2330cde88ab2c94c4e5b15915a608)) - clara.bayley
- fix doxstring errors - ([5a02ed9](https://github.com/yoctoyotta1024/CLEO/commit/5a02ed9cf9299c686b44dc7f673a9e1a39558b8e)) - clara.bayley
- note on installing mpi4py on levante - ([fa4d8f4](https://github.com/yoctoyotta1024/CLEO/commit/fa4d8f4fa4943bffaf1c11f58dc1a041ed69c355)) - clara.bayley
- correct mamba/conda - ([0911c34](https://github.com/yoctoyotta1024/CLEO/commit/0911c34f901d5c29f030292fa750c99c1b27e180)) - clara.bayley
#### Performance Improvements
- two commands in one - ([fd8ea03](https://github.com/yoctoyotta1024/CLEO/commit/fd8ea03e500a3b76bb3853812ef515af3de69c26)) - clara.bayley
#### Refactoring
- cleaner env creation - ([b3f96a1](https://github.com/yoctoyotta1024/CLEO/commit/b3f96a1567729dfafb926313b7e527ad6f1b8851)) - clara.bayley

- - -

## [v0.38.1](https://github.com/yoctoyotta1024/CLEO/compare/107f77c58aa94206f3ba90a515f8927aec8cb980..v0.38.1) - 2025-03-24
#### Bug Fixes
- update openmpi package for levante when usign intel compiler - ([107f77c](https://github.com/yoctoyotta1024/CLEO/commit/107f77c58aa94206f3ba90a515f8927aec8cb980)) - clara.bayley

- - -

## [v0.38.0](https://github.com/yoctoyotta1024/CLEO/compare/0fe36e7ee4d460a91b8ca8ae46b1cf6a11cda008..v0.38.0) - 2025-03-21
#### Features
- update version number to current v0.38.0 - ([383ae3c](https://github.com/yoctoyotta1024/CLEO/commit/383ae3c0a7bb4c7fbb637afff7fcbdfbe8b08131)) - clara.bayley
#### Performance Improvements
- remove redundant variable - ([445a6d9](https://github.com/yoctoyotta1024/CLEO/commit/445a6d9845026abc898244b38529504f83a9f0a5)) - clara.bayley
#### Refactoring
- add option to color 2d motion plot - ([02b4579](https://github.com/yoctoyotta1024/CLEO/commit/02b457955d3937e82250bfe19f90985ff52d409d)) - clara.bayley
- add option to plot 2-d motion of already chosen superdrops - ([db63f5b](https://github.com/yoctoyotta1024/CLEO/commit/db63f5bea94bb41b6987e3c9eafa4a0432e1373a)) - clara.bayley
- add better option to plot without smoothing - ([0c64f7a](https://github.com/yoctoyotta1024/CLEO/commit/0c64f7a38f4a88155c1b2712b6c2b7ad58e4d7e4)) - clara.bayley
- create parent directories if not already existing - ([4892b40](https://github.com/yoctoyotta1024/CLEO/commit/4892b40b923dc78b38db2f5fd80506cefe4cb9f3)) - clara.bayley
- add option to plot on given fig,ax - ([0fe36e7](https://github.com/yoctoyotta1024/CLEO/commit/0fe36e7ee4d460a91b8ca8ae46b1cf6a11cda008)) - clara.bayley

- - -

## [v0.37.0](https://github.com/yoctoyotta1024/CLEO/compare/0660ceba4f48b21d60d8c5c5efd69d0bce61672e..v0.37.0) - 2025-03-20
#### Documentation
- clearer docstrings about urbg random number ranges - ([2468457](https://github.com/yoctoyotta1024/CLEO/commit/2468457e6e44957891f684e994332b20ed7174b8)) - clara.bayley
- correct docstring - ([0660ceb](https://github.com/yoctoyotta1024/CLEO/commit/0660ceba4f48b21d60d8c5c5efd69d0bce61672e)) - clara.bayley
#### Features
- new file for shuffling superdroplets algorithm - ([82bd2a2](https://github.com/yoctoyotta1024/CLEO/commit/82bd2a254b91534e42911dc68418299cad032127)) - clara.bayley
#### Refactoring
- more uses KCS::team_size instead of Kokkos::AUTO - ([2a000c1](https://github.com/yoctoyotta1024/CLEO/commit/2a000c14a2a71b1c0c65342ea094306b40a0c5e5)) - clara.bayley
- split shuffle implementations into .cpp file - ([ad01b1b](https://github.com/yoctoyotta1024/CLEO/commit/ad01b1b36091d2a3d59dcadfdfa26575018c072a)) - clara.bayley

- - -

## [v0.36.2](https://github.com/yoctoyotta1024/CLEO/compare/9d56849e3964350f79e5b76cf6253043199e15f6..v0.36.2) - 2025-03-14
#### Bug Fixes
- add cap to ventilation factor for droplet radii > ~3mm - ([9d56849](https://github.com/yoctoyotta1024/CLEO/commit/9d56849e3964350f79e5b76cf6253043199e15f6)) - clara.bayley

- - -

## [v0.36.1](https://github.com/yoctoyotta1024/CLEO/compare/c01d9b1dae3c71ff9cfe77d93114c7750b303dd1..v0.36.1) - 2025-03-12
#### Bug Fixes
- mistake in range of valid values for fisher-yates shuffle - ([0ed7d10](https://github.com/yoctoyotta1024/CLEO/commit/0ed7d10bfd96614886f1f88a6923b51ed02874d3)) - clara.bayley
#### Performance Improvements
- delete unused header file from examples - ([c01d9b1](https://github.com/yoctoyotta1024/CLEO/commit/c01d9b1dae3c71ff9cfe77d93114c7750b303dd1)) - clara.bayley
#### Refactoring
- kokkos configuration parameters in a struct - ([dd52dd7](https://github.com/yoctoyotta1024/CLEO/commit/dd52dd7157df52e482f7a062ecf1e3916632d4b5)) - clara.bayley
- inclue ventillation coefficient in condensation.evaporation equation - ([d489a5a](https://github.com/yoctoyotta1024/CLEO/commit/d489a5a19bd244d6ece2fafa3e3b4b008eb237c8)) - clara.bayley

- - -

## [v0.36.0](https://github.com/yoctoyotta1024/CLEO/compare/e06bdbabbf5e0d96994d4391b1f67403a1f83bb0..v0.36.0) - 2025-03-04
#### Bug Fixes
- avoid casting -1, use constants instead - ([4b6b42e](https://github.com/yoctoyotta1024/CLEO/commit/4b6b42ea84bc49e2f410ed8fbc71e18bfbce803e)) - clara.bayley
- make privacy of functions and capture of values compatible with GPUs - ([fa95f75](https://github.com/yoctoyotta1024/CLEO/commit/fa95f7573375e88a0995f41ecc901b069ce2eb81)) - clara.bayley
#### Features
- new concept to define acceptable structures for boundary conditions in superdroplet movement - ([3a8860b](https://github.com/yoctoyotta1024/CLEO/commit/3a8860b97fa73b0c7fddb218554231e8b9934e4e)) - clara.bayley
- new concept to define acceptable structures for transporting superdroplets around the domain - ([feadc96](https://github.com/yoctoyotta1024/CLEO/commit/feadc9665be9e6cdf3281a8b58886945233a7a0b)) - clara.bayley
#### Refactoring
- update timing for speedtest example solution - ([9fb90a2](https://github.com/yoctoyotta1024/CLEO/commit/9fb90a2c133d6af5b2e7846c58bf0b210ae8c2f2)) - clara.bayley
- faster cuda implementation with class capture - ([12fd7e1](https://github.com/yoctoyotta1024/CLEO/commit/12fd7e1599f3f88abb5d75b010d9b4b335171c7e)) - clara.bayley
- move functors outside of DoCondensation for clearer encapsulation - ([52402d6](https://github.com/yoctoyotta1024/CLEO/commit/52402d691ef616c68936602d67de3fbcbcfaaf9c)) - clara.bayley
- move functors outside of DoCollisions for better encapsulation and clarity - ([202bdea](https://github.com/yoctoyotta1024/CLEO/commit/202bdeaad33f9b9ec71485f0913a8f1b0d714102)) - clara.bayley
- capture by value in single thread shuffling - ([ea3d054](https://github.com/yoctoyotta1024/CLEO/commit/ea3d05426f526c9f721d1ae4a4a97befa263e881)) - clara.bayley
- capture by value in lambda for runcleo references - ([611cc6b](https://github.com/yoctoyotta1024/CLEO/commit/611cc6b06cd8de7cc028cc5f616a288aba82776d)) - clara.bayley
- move functors outside of SDMMethods for clearer encapsulation - ([26fe16a](https://github.com/yoctoyotta1024/CLEO/commit/26fe16a4483710bb5705ed7cc91cd70d1e8f11b7)) - clara.bayley
- use concepts in runcleo and examples to constrain boundary conditions and transport templates - ([02c3026](https://github.com/yoctoyotta1024/CLEO/commit/02c3026cdf2575a501b0d7a1bee61980a4b9e307)) - clara.bayley
- move functors outside of MoveSupersInDomain for better encapsulation - ([91321a8](https://github.com/yoctoyotta1024/CLEO/commit/91321a8bde566bdf61c932a58c742ba89fe6c42b)) - clara.bayley
- don't use set refs with team_member when not needing nested loops - ([ae6ae29](https://github.com/yoctoyotta1024/CLEO/commit/ae6ae29c7e7daac1937b09f3ffe17db688f91067)) - clara.bayley
- don't hold subview returned as unused variable - ([14f6210](https://github.com/yoctoyotta1024/CLEO/commit/14f6210a9e01fbb09d283f5e62b2d81bf3d11bf0)) - clara.bayley
- don't assert total nullsupers from collisions - ([e06bdba](https://github.com/yoctoyotta1024/CLEO/commit/e06bdbabbf5e0d96994d4391b1f67403a1f83bb0)) - clara.bayley

- - -

## [v0.35.1](https://github.com/yoctoyotta1024/CLEO/compare/41b23c3e0d0b2b55a7dfac6f1c76963b20ad2568..v0.35.1) - 2025-03-03
#### Bug Fixes
- export yac and yaxt root to names YAC_ROOT and YAXT_ROOT - ([1670bf3](https://github.com/yoctoyotta1024/CLEO/commit/1670bf3ee9d1d325661e75c84971a0a328e2c20a)) - clara.bayley
- bash syntax errors - ([41b23c3](https://github.com/yoctoyotta1024/CLEO/commit/41b23c3e0d0b2b55a7dfac6f1c76963b20ad2568)) - clara.bayley
#### Continuous Integration
- add coupled_dynamics and domain flags to builds - ([c78e993](https://github.com/yoctoyotta1024/CLEO/commit/c78e9930e922c973836b35cd74f6de474d20f2c2)) - clara.bayley
#### Refactoring
- set flags in bash scripts for cleo cmake macros using CLEO_BUILD_FLAGS - ([25ecbbe](https://github.com/yoctoyotta1024/CLEO/commit/25ecbbe28b022359d687e1ec588ce07843e0678b)) - clara.bayley
- improve build status messages - ([28eed21](https://github.com/yoctoyotta1024/CLEO/commit/28eed21e044cd281a160683b9eee10b497372ddb)) - clara.bayley
- rename yac root flags and delete ENABLE_YAC_COUPLING redundant flag - ([ef669ca](https://github.com/yoctoyotta1024/CLEO/commit/ef669ca0104b0438d5b663c0aa665de1dedd9a9c)) - clara.bayley
- CLEO_[XXX] macros (flags) in CMakeLists.txt to not build exmaples and/or roughpaper - ([87cb1b4](https://github.com/yoctoyotta1024/CLEO/commit/87cb1b427e1ace2cd181512b7956043cee1aa89c)) - clara.bayley
- CLEO_[XXX] macros (flags) in CMakeLists.txt which determine coupled_dynamics and domain type - ([cd348ea](https://github.com/yoctoyotta1024/CLEO/commit/cd348eadebc579a30eab6baa458035a4e3c0a376)) - clara.bayley

- - -

## [v0.35.0](https://github.com/yoctoyotta1024/CLEO/compare/d5f98316f312a19b9be9fc5f754fb6786515271b..v0.35.0) - 2025-02-28
#### Bug Fixes
- fix missing yacyaxtroot arg to bash scripts - ([0ff84d2](https://github.com/yoctoyotta1024/CLEO/commit/0ff84d217c51086d7a012785e79ad5ad44ecc717)) - yoctoyotta1024
#### Documentation
- **(examples)** docstring typo fix - ([d5f9831](https://github.com/yoctoyotta1024/CLEO/commit/d5f98316f312a19b9be9fc5f754fb6786515271b)) - Sylwester Arabas
- better file explanation - ([dd9fe6d](https://github.com/yoctoyotta1024/CLEO/commit/dd9fe6d7082a67d05d360d85a460a7962bbc0b3f)) - yoctoyotta1024
#### Features
- new bash script to run fromfile example on juwels - ([383343a](https://github.com/yoctoyotta1024/CLEO/commit/383343ae626afea6a3238ac66a2d5ed93867746e)) - yoctoyotta1024
#### Miscellaneous Chores
- correct typo in docstring - ([8ad406e](https://github.com/yoctoyotta1024/CLEO/commit/8ad406edce2487b12bc80e624868cec9f0cdd2f2)) - yoctoyotta1024
- add TODOs - ([42ae76d](https://github.com/yoctoyotta1024/CLEO/commit/42ae76dd1b47a6f225070bfa9552b91ea41a23f6)) - yoctoyotta1024
#### Performance Improvements
- formatting - ([7a9dc1f](https://github.com/yoctoyotta1024/CLEO/commit/7a9dc1fe2194f2bc27cd9fc0cf371511779643af)) - clara.bayley
- delete unwanted comment - ([28bea10](https://github.com/yoctoyotta1024/CLEO/commit/28bea10cd24acf853d3af1ed3d6e3220a55a4966)) - yoctoyotta1024
#### Refactoring
- file rename - ([9610dbe](https://github.com/yoctoyotta1024/CLEO/commit/9610dbe1ea413285a0cf6d43e11493ac95c0a119)) - clara.bayley
- move levante bash scripts into levante - ([f6b0d3e](https://github.com/yoctoyotta1024/CLEO/commit/f6b0d3edcbdaf02537955f4b659b392f3cef68e6)) - clara.bayley
- remove macros for non-gpu superdrop functions - ([81dc834](https://github.com/yoctoyotta1024/CLEO/commit/81dc834e06f94d3bde5259b347568ae52c5de431)) - yoctoyotta1024

- - -

## [v0.34.0](https://github.com/yoctoyotta1024/CLEO/compare/ea35ec8cbcbc24a252621ab656e15cf98f584e26..v0.34.0) - 2025-02-24
#### Bug Fixes
- make examples compatible with refactored movesupersindomain struct - ([dd9b56d](https://github.com/yoctoyotta1024/CLEO/commit/dd9b56d91834436d78b6f7e1915f99e8693a2d0b)) - clara.bayley
- make examples compatible with refactored movesupersindomain struct - ([2e4de85](https://github.com/yoctoyotta1024/CLEO/commit/2e4de8517e0bec5c60fc0ffac463d78be6774098)) - yoctoyotta1024
- fix missing dynamic libraries at runtime - ([41f32c5](https://github.com/yoctoyotta1024/CLEO/commit/41f32c54bde7c996bf804f2b881c952edeb5fabe)) - yoctoyotta1024
- remove spack unload - ([05e9774](https://github.com/yoctoyotta1024/CLEO/commit/05e9774812f50a46045708d7bf65f687f951c04e)) - yoctoyotta1024
- path to juwels bash folder and juwels_packages - ([b7616f1](https://github.com/yoctoyotta1024/CLEO/commit/b7616f1a02be6b2e7b1d6bca7c3552b4875b6d61)) - yoctoyotta1024
#### Features
- new bash script to run divfree2d on juwels - ([296ebb7](https://github.com/yoctoyotta1024/CLEO/commit/296ebb79d90aa4f813fe0860c4ae6f66d663b0ad)) - yoctoyotta1024
- new bash scripts to run CLEO on JUWELS - ([3b4619f](https://github.com/yoctoyotta1024/CLEO/commit/3b4619fbc9f64b5e3faa67133f730db811cb9641)) - yoctoyotta1024
#### Performance Improvements
- linting files - ([38d1b2f](https://github.com/yoctoyotta1024/CLEO/commit/38d1b2fa7e146cb596f93317aecb333a8665085d)) - clara.bayley
#### Refactoring
- rename files and move into movement directory - ([c6cb78e](https://github.com/yoctoyotta1024/CLEO/commit/c6cb78e22b743a6b42b6a79874810d09f627a0f0)) - yoctoyotta1024
- movement of superdroplets across domain in seperate structure to MoveSupersInDomain - ([b67f5d2](https://github.com/yoctoyotta1024/CLEO/commit/b67f5d2ed608dfe48bf38a055f393927832ea041)) - yoctoyotta1024
- delete redundant file - ([78791dc](https://github.com/yoctoyotta1024/CLEO/commit/78791dc91a47751ffec81a7e77460a1a04e5bb6a)) - yoctoyotta1024
- rename motion -> sdmotion for clarity - ([5b80779](https://github.com/yoctoyotta1024/CLEO/commit/5b8077906594a97193c834a9cca0dc26c9711d12)) - yoctoyotta1024
- use ParaStationMPI not OpenMPI with gcc compiler - ([f0f5a98](https://github.com/yoctoyotta1024/CLEO/commit/f0f5a9880f0f69bb61d465ceb11aa70ba68d0681)) - yoctoyotta1024
- change runtime settings - ([62c5cf5](https://github.com/yoctoyotta1024/CLEO/commit/62c5cf58bba6e08c65d24fdff15f4e09dc6935b8)) - yoctoyotta1024
- change SLURM settings - ([89b428b](https://github.com/yoctoyotta1024/CLEO/commit/89b428baecdbe471810b11e08af96491372b6d2e)) - yoctoyotta1024
- lower cpu count for compiling - ([3e9e954](https://github.com/yoctoyotta1024/CLEO/commit/3e9e954c086c3ee6e2779ddccd3b6c4a75783aa2)) - yoctoyotta1024
- change default path to CLEO repo - ([6a84aa3](https://github.com/yoctoyotta1024/CLEO/commit/6a84aa34ef5b4e0c65838a7090c69aa2848bc876)) - yoctoyotta1024
- use juwels packages for intel compilers - ([f92bf06](https://github.com/yoctoyotta1024/CLEO/commit/f92bf06d073d2ba454e0ebe5045cf5efd2cd4033)) - yoctoyotta1024
- don't support YAC builds on JUWELS - ([4e9967e](https://github.com/yoctoyotta1024/CLEO/commit/4e9967ebc76756e0c291cc83adf8f7c150f8cf3e)) - yoctoyotta1024
- don't support CUDA builds on JUWELS - ([43da6bf](https://github.com/yoctoyotta1024/CLEO/commit/43da6bf770f4a875d22a2696505ad05963c3f3c5)) - yoctoyotta1024
- use juwels packages for gcc compilers - ([a32a451](https://github.com/yoctoyotta1024/CLEO/commit/a32a4519915ca38076be6e244d15b66b8be24b31)) - yoctoyotta1024
- set team size for heirarchal parallelism - ([a5c381e](https://github.com/yoctoyotta1024/CLEO/commit/a5c381e751b2de15ecfaf6b0ee06eebe49eb6a09)) - clara.bayley
- set number of host threads in fromfile example config - ([9d833c7](https://github.com/yoctoyotta1024/CLEO/commit/9d833c7ed961b75766ff750f6dcd08db738d7ba8)) - clara.bayley
- fromfile takes ntasks as argument - ([ea35ec8](https://github.com/yoctoyotta1024/CLEO/commit/ea35ec8cbcbc24a252621ab656e15cf98f584e26)) - clara.bayley

- - -

## [v0.33.1](https://github.com/yoctoyotta1024/CLEO/compare/2b624c97d680e468a97db0e1873347cf7f8e06ec..v0.33.1) - 2025-01-30
#### Bug Fixes
- gcc compiler error from taking address of rvalue - ([71a58db](https://github.com/yoctoyotta1024/CLEO/commit/71a58db95c403b9a69d430ad190a6bbee5645ca3)) - clara.bayley
#### Refactoring
- use argparse for fromfile args - ([05cc17c](https://github.com/yoctoyotta1024/CLEO/commit/05cc17c7dea9a56bd26884f4ddd6a2102ed603b1)) - clara.bayley
- move fromfile plotting to seperate script - ([e0e45b0](https://github.com/yoctoyotta1024/CLEO/commit/e0e45b0e3d40190d564c798a1cc9355a6257caee)) - clara.bayley
- add booleans to fromfile example run script - ([2b624c9](https://github.com/yoctoyotta1024/CLEO/commit/2b624c97d680e468a97db0e1873347cf7f8e06ec)) - clara.bayley

- - -

## [v0.33.0](https://github.com/yoctoyotta1024/CLEO/compare/8f483c9f2dafc32dcd93ca2afb66a31ef97a21a9..v0.33.0) - 2025-01-28
#### Bug Fixes
- make GPU compatible - ([14d0562](https://github.com/yoctoyotta1024/CLEO/commit/14d0562901d321a510a5c0d6fb51880acba3f341)) - clara.bayley
- make find_domainrefs GPU compatible - ([17bfc0d](https://github.com/yoctoyotta1024/CLEO/commit/17bfc0df4eea6aa8e60b24937c13034f088db17c)) - clara.bayley
#### Features
- add more plugs to profile superdroplet motion - ([accfbff](https://github.com/yoctoyotta1024/CLEO/commit/accfbff1eb3519d246a456cad45a098edee66151)) - clara.bayley
- use gbxs in create_cumlcounts function to avoid atomic conflicts - ([aa1c31f](https://github.com/yoctoyotta1024/CLEO/commit/aa1c31fe30ede49c1aa9c0c906b4bd787a4f7447)) - clara.bayley
- use gbxs in counting sort algorithm to reduce atomic conflicts - ([042e4f5](https://github.com/yoctoyotta1024/CLEO/commit/042e4f5fc0dc4a1c7cb6a86660f0b05c3be78849)) - clara.bayley
#### Refactoring
- use functor for create_cumlcounts loop - ([5140f0c](https://github.com/yoctoyotta1024/CLEO/commit/5140f0c18dd75e9864a342d003f7f59d3e24e34c)) - clara.bayley
- remove optional extras from sorting algorithm - ([3345c04](https://github.com/yoctoyotta1024/CLEO/commit/3345c0467f8b0e66db88a58e15014f65941ed8a5)) - clara.bayley
- use scatter view for counts summation to abstract atomics - ([efb47c5](https://github.com/yoctoyotta1024/CLEO/commit/efb47c577749b76bc95d72de1f043ddb01b63ac9)) - clara.bayley
- use find_partition_point also in find_ref for outer level parallelism cases - ([8f483c9](https://github.com/yoctoyotta1024/CLEO/commit/8f483c9f2dafc32dcd93ca2afb66a31ef97a21a9)) - clara.bayley

- - -

## [v0.32.0](https://github.com/yoctoyotta1024/CLEO/compare/e3164e2e0059fb518492332850c05abf7f20c832..v0.32.0) - 2025-01-23
#### Features
- replace kokkos/std sorting algorithm with counting sort algorithm - ([491f841](https://github.com/yoctoyotta1024/CLEO/commit/491f841fb8af867ffbbe159d3dcff28ead91b593)) - clara.bayley
#### Refactoring
- assume first position in totsupers is start of in domain supers and add docstrings - ([e3164e2](https://github.com/yoctoyotta1024/CLEO/commit/e3164e2e0059fb518492332850c05abf7f20c832)) - clara.bayley

- - -

## [v0.31.0](https://github.com/yoctoyotta1024/CLEO/compare/1a24d6319861ff6b74b25d9d86508a6a218b3a98..v0.31.0) - 2025-01-23
#### Bug Fixes
- encapsulation of supers in parallel regions - ([61ed09a](https://github.com/yoctoyotta1024/CLEO/commit/61ed09ad62ad1cc2b790ca05cfecf52b9cc38a90)) - clara.bayley
- typos in bash script - ([e61fc3d](https://github.com/yoctoyotta1024/CLEO/commit/e61fc3d0931e884e78e5ece055cfee5c81eec5e0)) - clara.bayley
- correct initial conditions and method to get size of supers view - ([b543045](https://github.com/yoctoyotta1024/CLEO/commit/b5430452c9dcc89c53114a00548a06293e422e36)) - clara.bayley
#### Features
- new struct to handle domain superdroplets - ([5c62ef1](https://github.com/yoctoyotta1024/CLEO/commit/5c62ef18040daeffc44dc4ce405164d12063a782)) - clara.bayley
#### Miscellaneous Chores
- rename SupersInDomain object - ([d2529be](https://github.com/yoctoyotta1024/CLEO/commit/d2529be2296d99948c84b047fe875ae4203a72a1)) - clara.bayley
- rename observers supers sub-view - ([26e2eba](https://github.com/yoctoyotta1024/CLEO/commit/26e2ebae00b38d437ff87a501583522519ea369a)) - clara.bayley
- use auto - ([7d0fd95](https://github.com/yoctoyotta1024/CLEO/commit/7d0fd954eeabfb3ef48f9dc88502694761cc9e12)) - clara.bayley
#### Performance Improvements
- nicer expression to reference gbx - ([67036d3](https://github.com/yoctoyotta1024/CLEO/commit/67036d3abe282b92fbd27b06b5308d1c9fa584a5)) - clara.bayley
#### Refactoring
- control sorting of totsupers from inside SupersInDomain struct - ([1ac934a](https://github.com/yoctoyotta1024/CLEO/commit/1ac934afa5fcc6a1b829b1ea39516e77f1f17d81)) - clara.bayley
- use kokkos style element access and add assert for supers size - ([68a1aed](https://github.com/yoctoyotta1024/CLEO/commit/68a1aed168e7d900021572549360f5671519e1df)) - clara.bayley
- return totsupers after sorting - ([944fb5b](https://github.com/yoctoyotta1024/CLEO/commit/944fb5bfe83f7589bfb81c8b72ddf1b54bfdbf5c)) - clara.bayley
- initconds use gbxmaps for nullgbxs - ([c5e0268](https://github.com/yoctoyotta1024/CLEO/commit/c5e0268e94de1fff10c19d7f7852cadbbc95364a)) - clara.bayley
- remove supers entire view from supersingbx object - ([d8ebeb5](https://github.com/yoctoyotta1024/CLEO/commit/d8ebeb51691a72327a2ee5ba1d966aadded836f6)) - clara.bayley
- use domainsupers to alter superdroplets during motion - ([8a9d3e9](https://github.com/yoctoyotta1024/CLEO/commit/8a9d3e91ef4c29cb7c9b8e22c51064109725ce40)) - clara.bayley
- don't use view in predcorr deltas - ([1a24d63](https://github.com/yoctoyotta1024/CLEO/commit/1a24d6319861ff6b74b25d9d86508a6a218b3a98)) - clara.bayley

- - -

## [v0.30.1](https://github.com/yoctoyotta1024/CLEO/compare/ffbd7900b9ef26749c0fb9c6b619584b83c143a6..v0.30.1) - 2024-12-21
#### Bug Fixes
- add spdtest results for new bash scripts - ([ffbd790](https://github.com/yoctoyotta1024/CLEO/commit/ffbd7900b9ef26749c0fb9c6b619584b83c143a6)) - clara.bayley

- - -

## [v0.30.0](https://github.com/yoctoyotta1024/CLEO/compare/1eb2f8afdf612f148fbde1ced64bef78a737f0ec..v0.30.0) - 2024-12-21
#### Bug Fixes
- correctly pass stacksize_limit - ([3a90e8f](https://github.com/yoctoyotta1024/CLEO/commit/3a90e8ff3a3625fc8894cde4c5474ed4f878e92f)) - clara.bayley
- fix partition for building gpu - ([43f0467](https://github.com/yoctoyotta1024/CLEO/commit/43f04678f3dd03d4c65aa9bb3e95cccdcaa337f1)) - clara.bayley
- debugging new scripts and tinkering - ([069e504](https://github.com/yoctoyotta1024/CLEO/commit/069e504813c9a6b00324681851cfe8a9cb5a5e2b)) - clara.bayley
#### Features
- add runtime optimisations for Levante - ([dc88aed](https://github.com/yoctoyotta1024/CLEO/commit/dc88aed3e0c15bc34673a373fbbff271a3d372a0)) - clara.bayley
#### Miscellaneous Chores
- move files - ([c9fa5d4](https://github.com/yoctoyotta1024/CLEO/commit/c9fa5d4698d8385b019fcb1656d1a40001eb10a0)) - clara.bayley
#### Performance Improvements
- be clearer on arg descriptions - ([2be2e3b](https://github.com/yoctoyotta1024/CLEO/commit/2be2e3bd38c5a1fd70910537fc8d2f3c21c8f9f8)) - clara.bayley
#### Refactoring
- add YAC runtime settings - ([7873891](https://github.com/yoctoyotta1024/CLEO/commit/7873891ddb92f2823d436e32163dbdf65dba0b65)) - clara.bayley
- run examples with intel compiler unless cuda build - ([a61c2dc](https://github.com/yoctoyotta1024/CLEO/commit/a61c2dc68ecc0c8d1a875436e173edd3219730de)) - clara.bayley
- add intel compiler option - ([91ec4d5](https://github.com/yoctoyotta1024/CLEO/commit/91ec4d56c0fc0f54d9059bc937ea809e5efcbd75)) - clara.bayley
- update gcc compiler version and flags - ([1253918](https://github.com/yoctoyotta1024/CLEO/commit/125391819c9eb2cfdce849283c58858a21068e6e)) - clara.bayley
- move packages into seperate file - ([767c7e2](https://github.com/yoctoyotta1024/CLEO/commit/767c7e2b1492b4447c68a6cb745a2663b67b4f91)) - clara.bayley
- allow examples and run_cleo script to use same runtime settings - ([6f0fc26](https://github.com/yoctoyotta1024/CLEO/commit/6f0fc26b843aa2592d9614617850d8661cebe753)) - clara.bayley
- update slurm of running example submission scripts - ([500d1d0](https://github.com/yoctoyotta1024/CLEO/commit/500d1d04bb6de8220d5d0784ab41bf79e60c037b)) - clara.bayley
- failed exit to running examples if wrong name used - ([482e8f9](https://github.com/yoctoyotta1024/CLEO/commit/482e8f976ca6b63fe7ecad84794b81fe6b098211)) - clara.bayley
- first draft import from check_inputs function script - ([d120920](https://github.com/yoctoyotta1024/CLEO/commit/d120920e7a32b0d72b88df2084aafaa6b157458f)) - clara.bayley
- new scripts for compiling and running cleo first draft - ([9b5277d](https://github.com/yoctoyotta1024/CLEO/commit/9b5277d94f2ab0e454faf805afc97eeea8196bb8)) - clara.bayley
- move install yac helper script - ([87e6a94](https://github.com/yoctoyotta1024/CLEO/commit/87e6a940ea09ecb89e9311e8b0c39d14ad89be2a)) - clara.bayley
- delete old redundant bash build helper files - ([e19b12d](https://github.com/yoctoyotta1024/CLEO/commit/e19b12d3abd146f491bd2cbc11a562e2b443cc1f)) - clara.bayley
- new scripts for building cleo firsrt draft - ([3129682](https://github.com/yoctoyotta1024/CLEO/commit/312968233f9c0fd9a71524480a0b5af265d932eb)) - clara.bayley
- new bash script for interface to building and compiling CLEO - ([079b0a0](https://github.com/yoctoyotta1024/CLEO/commit/079b0a0c4c3a9c7b3c570322dd0f516cf1c5586d)) - clara.bayley
- move install yac helper script - ([1a6c932](https://github.com/yoctoyotta1024/CLEO/commit/1a6c9328d7393097ec79989a41497a0726c4dec4)) - clara.bayley
- modify inputfiles slurm settings - ([14ea428](https://github.com/yoctoyotta1024/CLEO/commit/14ea428b4e7f72deea45ae71355667b0473051ab)) - clara.bayley
- initialise kokkos from struct given by config - ([270c5ed](https://github.com/yoctoyotta1024/CLEO/commit/270c5ed2b8925d3571da7617b2183a1499b5f365)) - clara.bayley
- specify resoures for examples - ([78925f9](https://github.com/yoctoyotta1024/CLEO/commit/78925f9db5b965eec6e1e9bfff7cf7006eba6ac0)) - clara.bayley
- initialise kokkos from struct given by config - ([1eb2f8a](https://github.com/yoctoyotta1024/CLEO/commit/1eb2f8afdf612f148fbde1ced64bef78a737f0ec)) - clara.bayley

- - -

## [v0.29.5](https://github.com/yoctoyotta1024/CLEO/compare/283e6feaf54436720704f52736f09a3b348e934a..v0.29.5) - 2024-12-19
#### Bug Fixes
- revert parallelising finding partition algorithm - ([c2b01e0](https://github.com/yoctoyotta1024/CLEO/commit/c2b01e0cec57546ff952ab8c0862fb79ae269dc4)) - clara.bayley
#### Documentation
- add note on experimental parallel version of find_partition_point - ([a6c126f](https://github.com/yoctoyotta1024/CLEO/commit/a6c126f2a956088021981a10586e57ff3dd186f4)) - clara.bayley
#### Miscellaneous Chores
- delete redundant functions - ([7df6112](https://github.com/yoctoyotta1024/CLEO/commit/7df6112bbb2af157d57b443a74a0d8cd56714278)) - clara.bayley
- add note on paths in bash script - ([4bc4dad](https://github.com/yoctoyotta1024/CLEO/commit/4bc4dad18fda3b8ddb2b2fc33755f27aff604517)) - clara.bayley
- use auto - ([75d9ad0](https://github.com/yoctoyotta1024/CLEO/commit/75d9ad0a3d19dd60f62cef5f514356d5d110f9d0)) - clara.bayley
- formatting and use auto - ([cd6f393](https://github.com/yoctoyotta1024/CLEO/commit/cd6f3938a50ec8450cefc0514a4fe045bb57f033)) - clara.bayley
- use auto - ([44dc3ae](https://github.com/yoctoyotta1024/CLEO/commit/44dc3ae35000715331e8267f5b0540bf47afbf8c)) - clara.bayley
- update kokkos version - ([283e6fe](https://github.com/yoctoyotta1024/CLEO/commit/283e6feaf54436720704f52736f09a3b348e934a)) - clara.bayley
#### Performance Improvements
- use kokkos min function not selfmade one - ([6e83288](https://github.com/yoctoyotta1024/CLEO/commit/6e83288d56fac859f4d43f4e3be2c07a57fd6c4f)) - clara.bayley
#### Refactoring
- new spdtest results for performance comparison - ([83998c0](https://github.com/yoctoyotta1024/CLEO/commit/83998c09ec372840aa0137c899ef4f8d45394c3a)) - clara.bayley
- rename gbxmaps ndims - ([00d527b](https://github.com/yoctoyotta1024/CLEO/commit/00d527ba0f2d306d281b107a065dae443ace828d)) - clara.bayley
- move setting of oob_gbxindex key in maps to optimised function - ([429d371](https://github.com/yoctoyotta1024/CLEO/commit/429d3713273b0079908615bd185c35899fc1054d)) - clara.bayley
- optimise null maps initialisation - ([3e108d7](https://github.com/yoctoyotta1024/CLEO/commit/3e108d7e777cd3ed40d6383414f8be819a420f33)) - clara.bayley
- optimise map initialisation of 3D model - ([a2ed085](https://github.com/yoctoyotta1024/CLEO/commit/a2ed085ac595236ddd3ab08fcc978311ca0f6b66)) - clara.bayley
- edit cartesian maps names and constructor and use auto - ([d5a2478](https://github.com/yoctoyotta1024/CLEO/commit/d5a24780995ded2a3c1b75bf077c01423fa19fd1)) - clara.bayley
- nested parallelisism for iscorrect function - ([58b10cb](https://github.com/yoctoyotta1024/CLEO/commit/58b10cb7cc627bb606e8f9f56a015de7350e7647)) - clara.bayley
- replace invalid argument with cassert and parallelise checking of gridboxes - ([54e5fd5](https://github.com/yoctoyotta1024/CLEO/commit/54e5fd59d12d12971873bc9150a4d396cf52c25e)) - clara.bayley
- parallelise finding partition point for refs - ([03e471f](https://github.com/yoctoyotta1024/CLEO/commit/03e471fa374e529f31fad90fba48ae12cc4da6f5)) - clara.bayley
- improve performance of nested parallelism partition point finding algorithm - ([0a1cc06](https://github.com/yoctoyotta1024/CLEO/commit/0a1cc06a60bea2740650904923f4d5f55b48fdd6)) - clara.bayley
- add bool to prevent default print statements - ([7e4c910](https://github.com/yoctoyotta1024/CLEO/commit/7e4c9101a2fbc45d41e2d4e2ed270ad23658ed70)) - clara.bayley

- - -

## [v0.29.4](https://github.com/yoctoyotta1024/CLEO/compare/e86a380d3339617c30d1eab992f322ffee50b8b8..v0.29.4) - 2024-12-13
#### Bug Fixes
- correct call signature for shima init conds - ([ca75dca](https://github.com/yoctoyotta1024/CLEO/commit/ca75dca3c86e8d5d6036e8e19abe06f23a674942)) - clara.bayley
#### Documentation
- update info on speed test example - ([42400f4](https://github.com/yoctoyotta1024/CLEO/commit/42400f486450303d08b40fd548683e43b8f503a6)) - clara.bayley
#### Refactoring
- add option for savelabel to gbx and thermo plots - ([c2b227a](https://github.com/yoctoyotta1024/CLEO/commit/c2b227a78acb14a092f7dae734646a9773f547c1)) - clara.bayley
- delete stats_filename parameter of model - ([f98dd14](https://github.com/yoctoyotta1024/CLEO/commit/f98dd14378d26d82e68de25637539eceb65bbbf9)) - clara.bayley
- remove redundant stats_filename from examples config files - ([e86a380](https://github.com/yoctoyotta1024/CLEO/commit/e86a380d3339617c30d1eab992f322ffee50b8b8)) - clara.bayley

- - -

## [v0.29.3](https://github.com/yoctoyotta1024/CLEO/compare/5564bde07ed6603c079c815f9d5cc79613f55f33..v0.29.3) - 2024-12-10
#### Bug Fixes
- remove unused class capture from lambda - ([6e99b63](https://github.com/yoctoyotta1024/CLEO/commit/6e99b633df8b45c5a94a408c6c7f6afe6777fdc7)) - clara.bayley
#### Refactoring
- add missing c++ standard lib include - ([7394a8d](https://github.com/yoctoyotta1024/CLEO/commit/7394a8ddddca144b0861ca010555f0018aab080d)) - clara.bayley
- change interface to xiprobdist calc and add new class to set minimum value of any distribution - ([5564bde](https://github.com/yoctoyotta1024/CLEO/commit/5564bde07ed6603c079c815f9d5cc79613f55f33)) - clara.bayley

- - -

## [v0.29.2](https://github.com/yoctoyotta1024/CLEO/compare/8e8ecbf76cde7d01baa95bc350a4a92a9e183446..v0.29.2) - 2024-12-06
#### Bug Fixes
- create superdroplets at domain top with xi>=1 - ([4a7be0d](https://github.com/yoctoyotta1024/CLEO/commit/4a7be0d5e236dcc203b973c1c38776d1d35ee205)) - clara.bayley
#### Refactoring
- add bool in SD creation to prevent un-physical superdroplets by default - ([8e8ecbf](https://github.com/yoctoyotta1024/CLEO/commit/8e8ecbf76cde7d01baa95bc350a4a92a9e183446)) - clara.bayley

- - -

## [v0.29.1](https://github.com/yoctoyotta1024/CLEO/compare/56a3b9a996f2654669794892d6ea83018ca2d6ca..v0.29.1) - 2024-12-06
#### Bug Fixes
- reverse order of longitude edge centers - ([f496c6e](https://github.com/yoctoyotta1024/CLEO/commit/f496c6efb1e7528c992bc205c2242099d00b338b)) - clara.bayley
- update python use for yac and yac_cadd_interp_stack_config_nnn call for latest yac version - ([d268b7c](https://github.com/yoctoyotta1024/CLEO/commit/d268b7ca6f0a40cc5428d993994f15111be5f5e8)) - clara.bayley
#### Miscellaneous Chores
- update yac and yaxt versions on CI - ([fcf576d](https://github.com/yoctoyotta1024/CLEO/commit/fcf576d888b3c3e7384698d004035e67dafe1574)) - clara.bayley
#### Refactoring
- adapt cleo domain for bubble to crop inside central portion of icon domain - ([b8e7de4](https://github.com/yoctoyotta1024/CLEO/commit/b8e7de4f31a11d0519b27b6a2beed8ab307831f4)) - clara.bayley
- adapt cleo domain for bubble to match entire icon domain - ([37958a0](https://github.com/yoctoyotta1024/CLEO/commit/37958a09c4aa739ca654a4311567bd85f0449428)) - clara.bayley
- update install yac bash script to make python bindings correctly - ([7ee4a84](https://github.com/yoctoyotta1024/CLEO/commit/7ee4a84e0ba564d3c3560938032bc0348785d2ca)) - clara.bayley
- update yacyaxt root dir - ([9c163be](https://github.com/yoctoyotta1024/CLEO/commit/9c163be2e43953bdd6535719158f0a40aa1d7afa)) - clara.bayley
- update expected solutions from spdtest - ([4e53e5b](https://github.com/yoctoyotta1024/CLEO/commit/4e53e5bbc6bea1a44949432d1debe05964a3ad2f)) - clara.bayley
- delete run_stats observer - ([56a3b9a](https://github.com/yoctoyotta1024/CLEO/commit/56a3b9a996f2654669794892d6ea83018ca2d6ca)) - clara.bayley

- - -

## [v0.29.0](https://github.com/yoctoyotta1024/CLEO/compare/f30212e726ca35fe15489fc9a81bfa9d498da0f0..v0.29.0) - 2024-12-05
#### Features
- add kokkos profiling hooks to measure computational performance - ([f30212e](https://github.com/yoctoyotta1024/CLEO/commit/f30212e726ca35fe15489fc9a81bfa9d498da0f0)) - clara.bayley

- - -

## [v0.28.4](https://github.com/yoctoyotta1024/CLEO/compare/aba51cd720e11577a53c3a7947295f10c3307b3f..v0.28.4) - 2024-11-27
#### Bug Fixes
- bubble run script uses python pathlib - ([aba51cd](https://github.com/yoctoyotta1024/CLEO/commit/aba51cd720e11577a53c3a7947295f10c3307b3f)) - clara.bayley

- - -

## [v0.28.3](https://github.com/yoctoyotta1024/CLEO/compare/fbad594c2925992ae0feade58d08a983c30accfd..v0.28.3) - 2024-11-22
#### Bug Fixes
- revert running fromfile example for longer - ([6761ffb](https://github.com/yoctoyotta1024/CLEO/commit/6761ffb942791409fe7e7a90a52f4698c8646dab)) - clara.bayley
#### Continuous Integration
- added script to the CI step for comparing parallel run results - ([5d710f0](https://github.com/yoctoyotta1024/CLEO/commit/5d710f02c3f0470f9cba60d5ccfbcd898989eb41)) - Wilton Jaciel Loch
- added parallelization execution test to verify that parallel execution is possible - ([fbad594](https://github.com/yoctoyotta1024/CLEO/commit/fbad594c2925992ae0feade58d08a983c30accfd)) - Wilton Jaciel Loch
#### Miscellaneous Chores
- Merge branch 'parallelization_ci' of https://github.com/wiltonloch/CLEO into parallelization_ci - ([aa494b8](https://github.com/yoctoyotta1024/CLEO/commit/aa494b8b447abe490f26b90a0a681555bd1d49ad)) - clara.bayley
#### Performance Improvements
- replaced global communication in superdrops exchange by p2p calls - ([ba4463a](https://github.com/yoctoyotta1024/CLEO/commit/ba4463a8a0eb7287c46f2263124395a83afd2a31)) - Wilton Jaciel Loch

- - -

## [v0.28.2](https://github.com/yoctoyotta1024/CLEO/compare/b1d57e5d52a245e6ca3904f047dce7c4f2541ac3..v0.28.2) - 2024-11-22
#### Bug Fixes
- correct tarball link - ([fbf3a28](https://github.com/yoctoyotta1024/CLEO/commit/fbf3a285c812601d886718719d1168ddbf182934)) - clara.bayley
- revert running fromfile example for longer - ([0618278](https://github.com/yoctoyotta1024/CLEO/commit/0618278dda3798581dc9e2632868455de97b6058)) - clara.bayley
#### Continuous Integration
- added script to the CI step for comparing parallel run results - ([74fb42a](https://github.com/yoctoyotta1024/CLEO/commit/74fb42a5dd4509b645d76f204a23feba7c706938)) - Wilton Jaciel Loch
- added parallelization execution test to verify that parallel execution is possible - ([b1d57e5](https://github.com/yoctoyotta1024/CLEO/commit/b1d57e5d52a245e6ca3904f047dce7c4f2541ac3)) - Wilton Jaciel Loch

- - -

## [v0.28.1](https://github.com/yoctoyotta1024/CLEO/compare/aab95ef25bdf5d7712b092ad7e3e5bd89014e770..v0.28.1) - 2024-11-22
#### Bug Fixes
- update levante mamba env path - ([aab95ef](https://github.com/yoctoyotta1024/CLEO/commit/aab95ef25bdf5d7712b092ad7e3e5bd89014e770)) - clara.bayley

- - -

## [v0.28.0](https://github.com/yoctoyotta1024/CLEO/compare/d783df060fed33171cfd5024e1c6c33e7b45c80e..v0.28.0) - 2024-11-21
#### Bug Fixes
- renamed function - ([5597588](https://github.com/yoctoyotta1024/CLEO/commit/55975888f1dd590afc11f4776087cc6ad544c511)) - clara.bayley
- fix order of dataset includes - ([643588a](https://github.com/yoctoyotta1024/CLEO/commit/643588a02ba8d88bd854934f17e9c212db7a65aa)) - clara.bayley
- deltas in predcorr gpu compatible and add kokkos macros to gbxmaps gpu functions - ([d03e2ec](https://github.com/yoctoyotta1024/CLEO/commit/d03e2ecee1d07f74d3ea8046479dfefc26ac5539)) - clara.bayley
- add MPI guards to fromfile_irreg example - ([8f02d55](https://github.com/yoctoyotta1024/CLEO/commit/8f02d5512e1a90be6a8cc4502d60bb7b57e78933)) - clara.bayley
- change executable names in CI - ([c930271](https://github.com/yoctoyotta1024/CLEO/commit/c930271279f08c588917914017106f53df985f3a)) - clara.bayley
- typo in comment - ([43e08fd](https://github.com/yoctoyotta1024/CLEO/commit/43e08fdbcf42b9f7892ab3f4bc8f20187832e266)) - clara.bayley
- fix sphinx dependencies after sphinx version 8 - ([92e396b](https://github.com/yoctoyotta1024/CLEO/commit/92e396ba0ec9345ec1b7b7485cb079b064d90a9a)) - clara.bayley
- Security vulnerability - ([d26a02c](https://github.com/yoctoyotta1024/CLEO/commit/d26a02c4091ef870102739820625a0f81a235319)) - clara.bayley
- added CartesianMaps instance to receive_dynamics in yac coupling dynamics - ([2370edb](https://github.com/yoctoyotta1024/CLEO/commit/2370edb69f8051b2595ddf80e7e353f249fb1f2d)) - Wilton Jaciel Loch
- added CartesianMaps instance to cvode and null coupling dynamics - ([95598de](https://github.com/yoctoyotta1024/CLEO/commit/95598dedcd48ea1f08ecc56fcb9e98f59bf1ef9f)) - Wilton Jaciel Loch
#### Continuous Integration
- updated cmake compiler flag to use mpi wrappers - ([97f9438](https://github.com/yoctoyotta1024/CLEO/commit/97f94389f48a8516609d2249a1f93bfb5a630083)) - Wilton Jaciel Loch
- minor changes for building with more restrict compiler rules - ([ee2636a](https://github.com/yoctoyotta1024/CLEO/commit/ee2636a86ba320d79501cf3e9ed41f95ec95ee31)) - Wilton Jaciel Loch
- updated cmake compiler flag to use mpi wrappers - ([36641c3](https://github.com/yoctoyotta1024/CLEO/commit/36641c32199cafd749bf3d2105b16c228c9b53ff)) - Wilton Jaciel Loch
#### Documentation
- add notes on coupling function calls - ([e3746d5](https://github.com/yoctoyotta1024/CLEO/commit/e3746d52416f7961f4388049543f350037b0663b)) - clara.bayley
- add openmpi compiler wrappers to requirements - ([ec3a531](https://github.com/yoctoyotta1024/CLEO/commit/ec3a531bc728b5c1082477a267e275f306fb786a)) - clara.bayley
- figures for memory layout - ([44ff4bb](https://github.com/yoctoyotta1024/CLEO/commit/44ff4bb90e098096e482dd17dbde4e02d346c01c)) - clara.bayley
- more intro on memory layout - ([a8da2a3](https://github.com/yoctoyotta1024/CLEO/commit/a8da2a337ad3699881713354af13177f402e228e)) - clara.bayley
- figure for timestepping - ([b1072f6](https://github.com/yoctoyotta1024/CLEO/commit/b1072f6e1d383bc49255ddbb2c5423db7a3dc976)) - clara.bayley
- more intro on timestepping - ([fa638e8](https://github.com/yoctoyotta1024/CLEO/commit/fa638e8acf88120c8aa2737b9f1ce3b1a8f426da)) - clara.bayley
- rearrange landing page - ([06f66c5](https://github.com/yoctoyotta1024/CLEO/commit/06f66c51483316361dea90c832f0712bb9df5d18)) - clara.bayley
- update build and executable names - ([69a71c3](https://github.com/yoctoyotta1024/CLEO/commit/69a71c35f56ddee11b74435de2d882c3f40f2c79)) - clara.bayley
- figures for memory layout - ([889f35a](https://github.com/yoctoyotta1024/CLEO/commit/889f35a077698f7c5f6c5234708a0961d1a8d37d)) - clara.bayley
- more intro on memory layout - ([66b63ce](https://github.com/yoctoyotta1024/CLEO/commit/66b63ce8aadf5aef0998d0436ccdaf4cc3ec507d)) - clara.bayley
- figure for timestepping - ([734536e](https://github.com/yoctoyotta1024/CLEO/commit/734536edb5a674d170218c5fde38270fa478c7f5)) - clara.bayley
- more intro on timestepping - ([1c2d27e](https://github.com/yoctoyotta1024/CLEO/commit/1c2d27e55fb80b81ece3252052a7e6bdd80de4f2)) - clara.bayley
- rearrange landing page - ([0b8ef33](https://github.com/yoctoyotta1024/CLEO/commit/0b8ef33dbabb86396247c4cf75188ef2fea27792)) - clara.bayley
#### Features
- new bash to submit slurm for all examples - ([c51c0d8](https://github.com/yoctoyotta1024/CLEO/commit/c51c0d859bd8f6052ae036b0c8488e1094bb98eb)) - clara.bayley
- new pysd module to help with creating and ploting initial condition binary files - ([780ee6e](https://github.com/yoctoyotta1024/CLEO/commit/780ee6ecf989da5c271fae4ef34900dcf212a123)) - clara.bayley
- add fromfile_irreg to build CI check - ([eed98b0](https://github.com/yoctoyotta1024/CLEO/commit/eed98b0e3348eaf511d2cea2f98ffaa249a781d0)) - clara.bayley
- add cmake target for formfile_irreg example - ([481e675](https://github.com/yoctoyotta1024/CLEO/commit/481e675bc5a196ae8e29470f5a4e536a63d2e38b)) - clara.bayley
- new example for irregular grid version of fromfile example (for MPI devlopment) - ([a32a6f8](https://github.com/yoctoyotta1024/CLEO/commit/a32a6f84c21481c4dba4613faf140070a4b5f800)) - clara.bayley
- added a collect_global_array implementation for long unsigned int type - ([77a43dc](https://github.com/yoctoyotta1024/CLEO/commit/77a43dc8a0677badd5302688baf39476c0282b51)) - Wilton Jaciel Loch
- initial mpi parallelization - ([d783df0](https://github.com/yoctoyotta1024/CLEO/commit/d783df060fed33171cfd5024e1c6c33e7b45c80e)) - Wilton Jaciel Loch
#### Miscellaneous Chores
- **(version)** v0.27.0 - ([be63380](https://github.com/yoctoyotta1024/CLEO/commit/be633804938bb910bcc0210547d6f69a169f84bb)) - yoctoyotta1024
- **(version)** v0.26.0 - ([d6d4c3d](https://github.com/yoctoyotta1024/CLEO/commit/d6d4c3df47c73f0b5075661f6f9f89ec7189f76f)) - yoctoyotta1024
- **(version)** v0.25.1 - ([b415082](https://github.com/yoctoyotta1024/CLEO/commit/b415082f410628d814f873f29ca207d80cf310ad)) - yoctoyotta1024
- delete redundant unused functions - ([91c397d](https://github.com/yoctoyotta1024/CLEO/commit/91c397d1436e1bb46a28ab761acd2e38c68e75c2)) - clara.bayley
- fix spelling mistakes - ([d8ec48b](https://github.com/yoctoyotta1024/CLEO/commit/d8ec48b671e998fd6375e0bec3367c7f3db4d970)) - clara.bayley
- formatting - ([660ffd0](https://github.com/yoctoyotta1024/CLEO/commit/660ffd0d7b6674a6a0529eecad29ff93ad95bb8d)) - clara.bayley
- formatting - ([e515e14](https://github.com/yoctoyotta1024/CLEO/commit/e515e14cb845d9b4b1a75bcb788451149d1d6063)) - clara.bayley
#### Refactoring
- add gbxmaps functions to avoid use of non-gpu compatible domain decomposition in single process builds - ([18b91ee](https://github.com/yoctoyotta1024/CLEO/commit/18b91eefba324535ed3b3048702e3630cabbcfe9)) - clara.bayley
- consistent use of out of bounds gbxindex value from constants - ([fb4617c](https://github.com/yoctoyotta1024/CLEO/commit/fb4617c0c7a8704a5bb2292f313916669a12b047)) - clara.bayley
- ensure gbxmaps returns correct types - ([7cbb859](https://github.com/yoctoyotta1024/CLEO/commit/7cbb8593390824fd11e28cfb9fc6dd8a767f2387)) - clara.bayley
- rework gridboxmaps - ([120109a](https://github.com/yoctoyotta1024/CLEO/commit/120109a81c90ea45a2992351ed869ae8f32e5c2e)) - clara.bayley
- rework predmotion - ([f427cea](https://github.com/yoctoyotta1024/CLEO/commit/f427cea60bf9b0011207c5f9a829cb90756a5cc6)) - clara.bayley
- move send/recv supers into seperate function with guard on comms > 1 - ([3bf7df7](https://github.com/yoctoyotta1024/CLEO/commit/3bf7df70993a7892f5622eab3b248716aeb8aafd)) - clara.bayley
- run fromfile example for longer and with 4 tasks - ([69a8382](https://github.com/yoctoyotta1024/CLEO/commit/69a83828be9100d2331e51274c4a35a37fffd5cd)) - clara.bayley
- helper functions for total_local_gridboxes and total_global_gridboxes - ([bd290d5](https://github.com/yoctoyotta1024/CLEO/commit/bd290d5152fa0b80baaa1c1000cece3e4e1e6c56)) - clara.bayley
- generalise couplingcomms for any gridbox maps - ([0754c4a](https://github.com/yoctoyotta1024/CLEO/commit/0754c4ad4aa6bb26e01b48ab4e54b377e3ad8e65)) - clara.bayley
- better use auto in examples - ([87f7398](https://github.com/yoctoyotta1024/CLEO/commit/87f7398046d205095a69dc9632aee4e5c3b35877)) - clara.bayley
- add to bash scripts for Levante the use openmpi compiler wrappers - ([b6d7e14](https://github.com/yoctoyotta1024/CLEO/commit/b6d7e14ee3446f005a7d09b997e51aaf29fdd062)) - clara.bayley
- added MPI capabilities to all roughpaper programs - ([ad28dc7](https://github.com/yoctoyotta1024/CLEO/commit/ad28dc730788c97cc9a27434801835e55888e13b)) - clara.bayley
- use new python module - ([dcfa7d2](https://github.com/yoctoyotta1024/CLEO/commit/dcfa7d2da008be34c91e832d7f0891721497ad3a)) - clara.bayley
- seperate steps in build CI - ([ccb73be](https://github.com/yoctoyotta1024/CLEO/commit/ccb73be05269990243f26a5e0f522f80485120e4)) - clara.bayley
- use pathlib for Paths properly - ([04b5e73](https://github.com/yoctoyotta1024/CLEO/commit/04b5e7313b05a4589fa8e5bde69b9b39f017f089)) - clara.bayley
- change levante account in bash scripts - ([52e6ace](https://github.com/yoctoyotta1024/CLEO/commit/52e6ace74b67fb23cf3c6999463ad2f07bc1a5d5)) - clara.bayley
- format figures - ([e2210ef](https://github.com/yoctoyotta1024/CLEO/commit/e2210efb92cff995b6002e4703d664dc2faa2fe0)) - clara.bayley
- add zXxXy dimensions in print statement - ([11dc5fe](https://github.com/yoctoyotta1024/CLEO/commit/11dc5fe56174a73098091081eac193687ca9345d)) - clara.bayley
- set irregular gbx boudndaries - ([78b07aa](https://github.com/yoctoyotta1024/CLEO/commit/78b07aa0d4dea1506cd30f0ccd30b828d549fa11)) - clara.bayley
- rename example fromfile -> fromfile_irreg - ([2f2550d](https://github.com/yoctoyotta1024/CLEO/commit/2f2550db56a205de2a15d89403589d43925a32f5)) - clara.bayley
- examples renaming to get rid of bad use of capital letters - ([91af179](https://github.com/yoctoyotta1024/CLEO/commit/91af179a42e223ddbd9ab6cfd1eaa8e494a83ab6)) - clara.bayley
- improve pre-commit hooks - ([f6d28be](https://github.com/yoctoyotta1024/CLEO/commit/f6d28be779f9b1d7adda960aa353700db344d0f8)) - clara.bayley
- added check to avoid sequential examples to be run with more than one MPI process - ([219c3c8](https://github.com/yoctoyotta1024/CLEO/commit/219c3c8c9d43e90687ff55afc7f69033c833050c)) - Wilton Jaciel Loch
- added MPI capabilities to all examples - ([48018c2](https://github.com/yoctoyotta1024/CLEO/commit/48018c22e90219a5246c07405ef1b04f8982b67b)) - Wilton Jaciel Loch
- added guard to test whether the sequential dataset has been included to allow sequential examples to run normally - ([e2e51d4](https://github.com/yoctoyotta1024/CLEO/commit/e2e51d4b9757f9607798ebc9cdbaa29f279cdd64)) - Wilton Jaciel Loch

- - -

## [v0.27.0](https://github.com/yoctoyotta1024/CLEO/compare/8fe37bd7e095b58363a798381730f706aaace06a..v0.27.0) - 2024-11-08
#### Features
- new bash to submit slurm for all examples - ([9443af2](https://github.com/yoctoyotta1024/CLEO/commit/9443af2827dd49936fdd20f224a648f51e823d2b)) - clara.bayley
- new pysd module to help with creating and ploting initial condition binary files - ([5f3949c](https://github.com/yoctoyotta1024/CLEO/commit/5f3949cc5666a6c7d8a65e92ffa5d7be66716ceb)) - clara.bayley
#### Miscellaneous Chores
- formatting - ([49959bd](https://github.com/yoctoyotta1024/CLEO/commit/49959bdbf33654d2b0d79cb9bf94d14f1a82ed72)) - clara.bayley
#### Refactoring
- use new python module - ([818bd3b](https://github.com/yoctoyotta1024/CLEO/commit/818bd3b332cc4f593258671eb0a6e8f14d5862ce)) - clara.bayley
- seperate steps in build CI - ([198df4d](https://github.com/yoctoyotta1024/CLEO/commit/198df4d81f1a5a0c30b8e50e74c6d60b8db165c9)) - clara.bayley
- use pathlib for Paths properly - ([83bc5dd](https://github.com/yoctoyotta1024/CLEO/commit/83bc5dd57e273da159abfc940f64b491331ad542)) - clara.bayley
- change levante account in bash scripts - ([8fe37bd](https://github.com/yoctoyotta1024/CLEO/commit/8fe37bd7e095b58363a798381730f706aaace06a)) - clara.bayley

- - -

## [v0.26.0](https://github.com/yoctoyotta1024/CLEO/compare/d720d22535e3698bf5066e31f4d41c4d53cc3faa..v0.26.0) - 2024-09-11
#### Bug Fixes
- change executable names in CI - ([642b5e1](https://github.com/yoctoyotta1024/CLEO/commit/642b5e1295447dd610c016b335b1a1169e19f458)) - clara.bayley
- typo in comment - ([80ecf65](https://github.com/yoctoyotta1024/CLEO/commit/80ecf653d00fe37187461d60d4a3570427f68b6f)) - clara.bayley
#### Documentation
- update build and executable names - ([e7d358e](https://github.com/yoctoyotta1024/CLEO/commit/e7d358eb11372aca3dc241c3efd16f436024c64c)) - clara.bayley
#### Features
- add fromfile_irreg to build CI check - ([a684930](https://github.com/yoctoyotta1024/CLEO/commit/a684930c01679b1f82fe298d1cea38cf9d7c0ad3)) - clara.bayley
- add cmake target for formfile_irreg example - ([34786e2](https://github.com/yoctoyotta1024/CLEO/commit/34786e2be179de5c25c33b9be1fc4a5cc1fd6020)) - clara.bayley
- new example for irregular grid version of fromfile example (for MPI devlopment) - ([3b05380](https://github.com/yoctoyotta1024/CLEO/commit/3b05380b107307c9f371da408fd8ff76359d39cc)) - clara.bayley
#### Refactoring
- format figures - ([477b29c](https://github.com/yoctoyotta1024/CLEO/commit/477b29c055ad820802041d93cee1e6c1bc244e49)) - clara.bayley
- add zXxXy dimensions in print statement - ([3fe2737](https://github.com/yoctoyotta1024/CLEO/commit/3fe273762825810848242e6382db7c1261aacfd0)) - clara.bayley
- set irregular gbx boudndaries - ([e287cab](https://github.com/yoctoyotta1024/CLEO/commit/e287cabc2b113f84a2f9649b12c69d1889e6d57b)) - clara.bayley
- rename example fromfile -> fromfile_irreg - ([de332ae](https://github.com/yoctoyotta1024/CLEO/commit/de332ae97bed58d9b7a3d2fa3317161a9eb29612)) - clara.bayley
- examples renaming to get rid of bad use of capital letters - ([d720d22](https://github.com/yoctoyotta1024/CLEO/commit/d720d22535e3698bf5066e31f4d41c4d53cc3faa)) - clara.bayley

- - -

## [v0.25.1](https://github.com/yoctoyotta1024/CLEO/compare/8cf9fc790d9de7578dd2f2c9bd59bfeeaa3fae7b..v0.25.1) - 2024-09-04
#### Bug Fixes
- fix sphinx dependencies after sphinx version 8 - ([54bb0fa](https://github.com/yoctoyotta1024/CLEO/commit/54bb0fa76172579ab0140f084c94aba1109ff4cf)) - clara.bayley
- Security vulnerability - ([9ab37dc](https://github.com/yoctoyotta1024/CLEO/commit/9ab37dc81f6478a180159fc53e280a7a6b16cf23)) - clara.bayley
#### Miscellaneous Chores
- formatting - ([b6210eb](https://github.com/yoctoyotta1024/CLEO/commit/b6210eb2904c7c3f66597e8a7ad3f0299cd01f46)) - clara.bayley
#### Refactoring
- improve pre-commit hooks - ([4a854c3](https://github.com/yoctoyotta1024/CLEO/commit/4a854c36ccae0e463ad202238109ccf664c0a6e9)) - clara.bayley
- move yac_raw_data_to_target_array into receive_yac_field function - ([963774d](https://github.com/yoctoyotta1024/CLEO/commit/963774d02d62a31d19a82db4a9205837069ce900)) - clara.bayley
- avoid unneccesary use of named variables - ([8cf9fc7](https://github.com/yoctoyotta1024/CLEO/commit/8cf9fc790d9de7578dd2f2c9bd59bfeeaa3fae7b)) - clara.bayley

- - -

## [v0.25.0](https://github.com/yoctoyotta1024/CLEO/compare/55cd54fed2074aa32c44d8398d2437c24e9c8c84..v0.25.0) - 2024-08-18
#### Bug Fixes
- correctly slice icon data to account for vertical levels being top down not bottom up - ([981f486](https://github.com/yoctoyotta1024/CLEO/commit/981f486106c9411a8811b731e84dda75f1c92bb5)) - clara.bayley
- correct ordering of yac_raw_data into target_array for each variable - ([3a111f5](https://github.com/yoctoyotta1024/CLEO/commit/3a111f586417aa28bdbde3588a99d0cefbcc34dc)) - clara.bayley
- add missin types - ([f598954](https://github.com/yoctoyotta1024/CLEO/commit/f5989545e4354607891629f03baff29c510dbea6)) - clara.bayley
- minor variable renaming - ([67ca555](https://github.com/yoctoyotta1024/CLEO/commit/67ca5551a047c57abba91d2645f72618e993546d)) - clara.bayley
- ensure gbxidxs for plotting are ints - ([2e4a1a9](https://github.com/yoctoyotta1024/CLEO/commit/2e4a1a9324c23e83e7cc4e38ca26f09325184900)) - clara.bayley
- ensure gbxidxs for plotting are ints - ([1ff8812](https://github.com/yoctoyotta1024/CLEO/commit/1ff88120d5cdd98f51631d3d6344e39c1e5477cc)) - clara.bayley
- use correct grid in tmp run script - ([63e5fe2](https://github.com/yoctoyotta1024/CLEO/commit/63e5fe2fd45e74c88a94027eb173dec76b7ede5c)) - clara.bayley
- minor bug fixes for typos and naming types for args into script - ([70914cf](https://github.com/yoctoyotta1024/CLEO/commit/70914cf5e55d68c4198f4ec2fc9248815ae206fa)) - clara.bayley
- don't activate cleoenv during compilation - ([b2f2946](https://github.com/yoctoyotta1024/CLEO/commit/b2f2946ef03fd332d93a81d861375f217aad0a24)) - clara.bayley
- correct gridfile for bubble - ([dd09574](https://github.com/yoctoyotta1024/CLEO/commit/dd09574b47b6325ffcd3f51a5f302dd3d0fc20a8)) - clara.bayley
- correct python path for yac - ([55cd54f](https://github.com/yoctoyotta1024/CLEO/commit/55cd54fed2074aa32c44d8398d2437c24e9c8c84)) - clara.bayley
#### Features
- chose file  for ICONN data fom arguments passed to python script - ([259ab37](https://github.com/yoctoyotta1024/CLEO/commit/259ab37fc879921a83062df11959b71da8456759)) - clara.bayley
#### Miscellaneous Chores
- extend bubble example run time - ([33c9cbb](https://github.com/yoctoyotta1024/CLEO/commit/33c9cbb49d46da54c6fc466e2c4f8c97774e8003)) - clara.bayley
- formatting - ([bdc10d1](https://github.com/yoctoyotta1024/CLEO/commit/bdc10d14a106ffb4dea3556a59723a2c218e19ff)) - clara.bayley
- formatting - ([e0c2cae](https://github.com/yoctoyotta1024/CLEO/commit/e0c2cae69e565735898efea7ab4d19340b5a2dcc)) - clara.bayley
- formatting - ([33c97c9](https://github.com/yoctoyotta1024/CLEO/commit/33c97c9c852511e64827ac805db10718e7f5d535)) - clara.bayley
- formatting - ([e1df3b3](https://github.com/yoctoyotta1024/CLEO/commit/e1df3b324dbd203a8acaf9e51db376da686abd04)) - clara.bayley
#### Refactoring
- bubble with higher resolution settings - ([d5a026c](https://github.com/yoctoyotta1024/CLEO/commit/d5a026c228b1cbbb3ae578b33cc7f5192b3dcfde)) - clara.bayley
- set good params for bubble example - ([df9f137](https://github.com/yoctoyotta1024/CLEO/commit/df9f1377bc7f66626e10c655a51e98bea1161b47)) - clara.bayley
- reorganise functionality to remove case switch - ([e74a3f2](https://github.com/yoctoyotta1024/CLEO/commit/e74a3f237fd8411d318e4309957bdcd957df0eb7)) - clara.bayley
- demand that yac is given max and min longitude and latitude in config params - ([3e67b02](https://github.com/yoctoyotta1024/CLEO/commit/3e67b02ebebdcaa38a8a47c1b7ff2501764e910c)) - clara.bayley
- change grid for cleo in bubble3d example - ([f1d9d9f](https://github.com/yoctoyotta1024/CLEO/commit/f1d9d9f7a5d884229d34927305aa42d50d686380)) - clara.bayley
- reorganise tmp run script for clarity - ([e4a2bf7](https://github.com/yoctoyotta1024/CLEO/commit/e4a2bf72c19d4459f2f22ba5182a92cdf57ca7cf)) - clara.bayley
- make ICON grid name an input variable and aadded print statement - ([cbc9c21](https://github.com/yoctoyotta1024/CLEO/commit/cbc9c212ef8d1dc4a3b5265732b0e3f11881571c)) - clara.bayley
- make ICON grid name an input variable and aadded print statement - ([7ffa41d](https://github.com/yoctoyotta1024/CLEO/commit/7ffa41d07d6712d6fa4cb49476f815b144b6cb39)) - clara.bayley
- make yac_bubble_data_reader python script more general + formatting - ([8f8b6bf](https://github.com/yoctoyotta1024/CLEO/commit/8f8b6bf05ec487670a47d43cc9c71a8298239395)) - clara.bayley
- remove make clean from tmp file - ([dda8571](https://github.com/yoctoyotta1024/CLEO/commit/dda857127c289a240e291f2c19870d3e2102c1d1)) - clara.bayley
- change obsstep, coupling step and add another superdroplet observer - ([cd7fc9d](https://github.com/yoctoyotta1024/CLEO/commit/cd7fc9dbea79459ac45a25058a2d5664de1cbd8b)) - clara.bayley
- make yac script compatible with Torus Triangles gridfile - ([ad01628](https://github.com/yoctoyotta1024/CLEO/commit/ad0162814eecf1c25f0b646762c3d9add23cfdb9)) - clara.bayley

- - -

## [v0.24.0](https://github.com/yoctoyotta1024/CLEO/compare/aec716f2a6071ccd05e18c2406c3024d0e092932..v0.24.0) - 2024-08-15
#### Bug Fixes
- delete yac_version after moving python bindings into yac folder - ([a5e4bc7](https://github.com/yoctoyotta1024/CLEO/commit/a5e4bc7ef24fe8aeb05f185f6556e87ad2669dd4)) - clara.bayley
- debugging build ci fail with latest sphinx - ([3a7ea19](https://github.com/yoctoyotta1024/CLEO/commit/3a7ea19816cf672c1d3123eb6ce258bbbac75b87)) - clara.bayley
- debugging build ci fail with latest sphinx - ([be7706e](https://github.com/yoctoyotta1024/CLEO/commit/be7706eb10885b73e44b3627d6501d0e0a5a99a8)) - clara.bayley
- debugging build ci fail with latest sphinx - ([b2f9dd5](https://github.com/yoctoyotta1024/CLEO/commit/b2f9dd5baf5a51993d652c72143f4c39fe0c4784)) - clara.bayley
- docs fix typo - ([7303dc3](https://github.com/yoctoyotta1024/CLEO/commit/7303dc37688079bf0c30032935dd1fb56fd6e8a7)) - clara.bayley
- bubble tmp bash script working - ([402fe2a](https://github.com/yoctoyotta1024/CLEO/commit/402fe2aedab6857cf8f8f92f39c4cfff09d48db2)) - clara.bayley
- add note on how to compile without loading cleoenv - ([430bcd5](https://github.com/yoctoyotta1024/CLEO/commit/430bcd5f1921ce8b2632c800bb9886ea8722f64a)) - clara.bayley
#### Features
- raw bubble script working - ([53474be](https://github.com/yoctoyotta1024/CLEO/commit/53474bedd202ade5360de0a1fc8c15224ecbb625)) - clara.bayley
- new explicit bash for debugging bubble run - ([9ec7225](https://github.com/yoctoyotta1024/CLEO/commit/9ec722561d06c86b15aab958dd1337e3e0576fd9)) - clara.bayley
- add module purge into tmp bash - ([782c997](https://github.com/yoctoyotta1024/CLEO/commit/782c9970c90683cc5cb382b8b5b7f76da27a6fd1)) - clara.bayley
#### Miscellaneous Chores
- delet coments - ([d6382e6](https://github.com/yoctoyotta1024/CLEO/commit/d6382e6fe298dfa92aaabcb24caab1fe0745a967)) - clara.bayley
- delete unwanted comment - ([7657f23](https://github.com/yoctoyotta1024/CLEO/commit/7657f2359542dff96502741205a5fd21abdc808b)) - clara.bayley
- uncomment yaxt install - ([3ec8d71](https://github.com/yoctoyotta1024/CLEO/commit/3ec8d713362e4d3b8270858230e2fbcc23d58206)) - clara.bayley
- comment with example of yac install call - ([f872ed7](https://github.com/yoctoyotta1024/CLEO/commit/f872ed7b26b3de48f29a0febfc64af2c099df8c6)) - clara.bayley
- comment for includes needed to run example with yac - ([58dbab6](https://github.com/yoctoyotta1024/CLEO/commit/58dbab6b30cc5dfa0d491f050b3da8a32fda346e)) - clara.bayley
#### Refactoring
- remove verbose bash script in install yac - ([5c2445d](https://github.com/yoctoyotta1024/CLEO/commit/5c2445d645b5058fb8fceaee2311e77c2d7d2f99)) - clara.bayley
- generalise temporary bubble test script - ([d818018](https://github.com/yoctoyotta1024/CLEO/commit/d8180189cca27ccd27a19922eb21446a5e21a73b)) - clara.bayley
- all stages of bubble3d in tmp script - ([31fa121](https://github.com/yoctoyotta1024/CLEO/commit/31fa121eea81312edf84608f139354547ec01dd7)) - clara.bayley
- change source paths for examples - ([bb7df5c](https://github.com/yoctoyotta1024/CLEO/commit/bb7df5ce706f9dfa593c5fb3f468027aa22f384b)) - clara.bayley
- delete unwanted compiler flags for yac installation - ([e1d6efe](https://github.com/yoctoyotta1024/CLEO/commit/e1d6efeb105f0a5be1320ae70ba1a38c2ea4075c)) - clara.bayley
- use set in bash script to print commands - ([2b05f9d](https://github.com/yoctoyotta1024/CLEO/commit/2b05f9d2a2fca00ec24cf7c797b089b29b6b85d5)) - clara.bayley
- enable python bindings in yac build - ([ae26dee](https://github.com/yoctoyotta1024/CLEO/commit/ae26deebfe743152f53f1a6bb4509d54d2e73eab)) - clara.bayley
- update yac version and pythonpath export - ([aec716f](https://github.com/yoctoyotta1024/CLEO/commit/aec716f2a6071ccd05e18c2406c3024d0e092932)) - clara.bayley

- - -

## [v0.23.2](https://github.com/yoctoyotta1024/CLEO/compare/597b43cbd9bde146c122d8d94ae65b3c3afa9976..v0.23.2) - 2024-07-24
#### Bug Fixes
- corrected templating over store of monitoring mass moment xarrays - ([597b43c](https://github.com/yoctoyotta1024/CLEO/commit/597b43cbd9bde146c122d8d94ae65b3c3afa9976)) - clara.bayley

- - -

## [v0.23.1](https://github.com/yoctoyotta1024/CLEO/compare/605b536b4fa564568a92e76cd9cf1915689426ad..v0.23.1) - 2024-07-19
#### Bug Fixes
- remove units of totmass_condensed - ([605b536](https://github.com/yoctoyotta1024/CLEO/commit/605b536b4fa564568a92e76cd9cf1915689426ad)) - clara.bayley

- - -

## [v0.23.0](https://github.com/yoctoyotta1024/CLEO/compare/e9c9dad54d10db958df7f024044feef653426644..v0.23.0) - 2024-07-19
#### Features
- new monitor to observe massmoments of raindrop size distribution during motion and microphysics - ([23db516](https://github.com/yoctoyotta1024/CLEO/commit/23db5164bdb2126587b95d6f7692a3c5259ab42b)) - clara.bayley
- new monitors xarrays for writing rain mass moments to dataset - ([0eb00cb](https://github.com/yoctoyotta1024/CLEO/commit/0eb00cbc2aef629f21f659b8d7c7f30e16347348)) - clara.bayley
- new monitors views for calculating rain mass moments - ([e9c9dad](https://github.com/yoctoyotta1024/CLEO/commit/e9c9dad54d10db958df7f024044feef653426644)) - clara.bayley
#### Refactoring
- use new observer in roughpaper src main for example - ([7f587cb](https://github.com/yoctoyotta1024/CLEO/commit/7f587cbb93a7d1aca6c2f2f9fec568fed5b4841f)) - clara.bayley
- template over xarray type for monitoring mass moments - ([01dfd13](https://github.com/yoctoyotta1024/CLEO/commit/01dfd139f1da9b22e0634192d03b32837fd3c7f2)) - clara.bayley

- - -

## [v0.22.0](https://github.com/yoctoyotta1024/CLEO/compare/e80b45263ba73cbf0f933714d45bea13d35dd487..v0.22.0) - 2024-07-09
#### Continuous Integration
- removing now inexistent yac examples - ([be38da0](https://github.com/yoctoyotta1024/CLEO/commit/be38da091294569f5ee2e2810f8be926a4576b76)) - wiltonloch
#### Features
- files created for running bubble3d example through python scripts like for other examples - ([15aff21](https://github.com/yoctoyotta1024/CLEO/commit/15aff211e90112d95de69b71b5f43689f0314a32)) - clara.bayley
#### Miscellaneous Chores
- restore fromfile example's CMakeLists.txt - ([96e3dd0](https://github.com/yoctoyotta1024/CLEO/commit/96e3dd0f5666d98c08767a71a0feaf4d1925f3b9)) - clara.bayley
- restore fromfile example - ([0cd1722](https://github.com/yoctoyotta1024/CLEO/commit/0cd17222e8c85b2236056042f5f8e8c15351b7bc)) - clara.bayley
- update todo note - ([9ebee4d](https://github.com/yoctoyotta1024/CLEO/commit/9ebee4d45edac8c916016727d8bbf595725c4a8e)) - clara.bayley
- rename file for consistency - ([2fe57d6](https://github.com/yoctoyotta1024/CLEO/commit/2fe57d65a94e92baaed9a6b391a4aa9564613862)) - clara.bayley
- delete unused files - ([ebe4583](https://github.com/yoctoyotta1024/CLEO/commit/ebe4583ee8faff13ff966a8a491810daa24e3a83)) - clara.bayley
- format header - ([7f33833](https://github.com/yoctoyotta1024/CLEO/commit/7f33833e7f0df6173d0c380e7c78ca1797c204df)) - clara.bayley
- format header - ([f566239](https://github.com/yoctoyotta1024/CLEO/commit/f566239aaa6a00d4179d74383ae313b80155f38b)) - clara.bayley
- delete already deleted subdir from cmake - ([cdd81bc](https://github.com/yoctoyotta1024/CLEO/commit/cdd81bc70d681485611459a6a25229c74982f7d7)) - clara.bayley
- rename executable to be consistent with other examples - ([7d10046](https://github.com/yoctoyotta1024/CLEO/commit/7d1004679fb419da7b4cc505d7a0c2cd8bc53a00)) - clara.bayley
- rename bubble_3d to bubble3d to be consistent with other examples - ([026e9d2](https://github.com/yoctoyotta1024/CLEO/commit/026e9d2e2cec7682a3d6797828b6d7ef0c8aedf5)) - clara.bayley
- simplified yac example folder and renamed it to bubble_3d - ([8dc0692](https://github.com/yoctoyotta1024/CLEO/commit/8dc069248316a8ea3a727d25c0cce8a6a11055e6)) - wiltonloch
- removed all yac mentions in the fromfile example - ([8677e3c](https://github.com/yoctoyotta1024/CLEO/commit/8677e3c565d27bf5022a95b0e4f971f3edfc8371)) - wiltonloch
- moved fromfile example from inside yac to dedicated example folder - ([5198220](https://github.com/yoctoyotta1024/CLEO/commit/5198220f1b47ebf747490d1dd9d458532ef9828a)) - wiltonloch
- removed yac divfreemotion stale example - ([e80b452](https://github.com/yoctoyotta1024/CLEO/commit/e80b45263ba73cbf0f933714d45bea13d35dd487)) - wiltonloch
#### Refactoring
- add fromfile to CI builds - ([0fb1168](https://github.com/yoctoyotta1024/CLEO/commit/0fb11680293af85a5b91a32f17b5c70326230fba)) - clara.bayley
- delete unwanted example of fromfile 3d dynamics - ([ea9176d](https://github.com/yoctoyotta1024/CLEO/commit/ea9176dd40f75648bfec941d47d49356f3927bb0)) - clara.bayley

- - -

## [v0.21.0](https://github.com/yoctoyotta1024/CLEO/compare/d8908a94f520210b94193916644dceb2d7c1ac76..v0.21.0) - 2024-06-27
#### Bug Fixes
- fixed coupling dimension order and added enum labels for indexing - ([1dcb012](https://github.com/yoctoyotta1024/CLEO/commit/1dcb012192b3931db1da5b9dca9fdffcb520d1e9)) - wiltonloch
- added conversion factor to receive_yac_field to get data in the correct units for CLEO - ([e786de4](https://github.com/yoctoyotta1024/CLEO/commit/e786de4f762f0f08b926a9e08ff2cd461e5fc869)) - wiltonloch
- added explicit cpp std endl to avoid output srun's multicomponent output suppression - ([16a87b8](https://github.com/yoctoyotta1024/CLEO/commit/16a87b83f34eb7006bf0fe92e20935bb12be9255)) - wiltonloch
- added more robust comparison for NaN values in the yac coupling - ([b6feaf5](https://github.com/yoctoyotta1024/CLEO/commit/b6feaf502eb7b7426b02d4364e549d55298d906b)) - wiltonloch
#### Features
- **(config)** allowing partial boundaries to be defined for the CLEO YAC grid - ([d2bdcb7](https://github.com/yoctoyotta1024/CLEO/commit/d2bdcb748bc91650a63fbf9cc3c3043c1cdb0872)) - wiltonloch
#### Miscellaneous Chores
- use cleo project for binary and source directories in cmake - ([33e0a58](https://github.com/yoctoyotta1024/CLEO/commit/33e0a58dc915aa97f2b163e2c97c2fb404a67e80)) - clara.bayley
- tidy up includes - ([2ec727b](https://github.com/yoctoyotta1024/CLEO/commit/2ec727bc2df4b3617a4fd83a54ccc36c9e7d020c)) - clara.bayley
- correct includes and formatting - ([d8908a9](https://github.com/yoctoyotta1024/CLEO/commit/d8908a94f520210b94193916644dceb2d7c1ac76)) - clara.bayley
#### Refactoring
- renamed subroutine that receives all the fields from YAC - ([976fe21](https://github.com/yoctoyotta1024/CLEO/commit/976fe218dae19267cf806bf4fec85817950813d1)) - wiltonloch
- changed field and coupling timesteps to match the ones in ICON - ([b1fb996](https://github.com/yoctoyotta1024/CLEO/commit/b1fb996be81fbda84db5de9bf21b1e3c82a99a6b)) - wiltonloch
- changed python reader script component and grid names to match ICON's - ([eeb3475](https://github.com/yoctoyotta1024/CLEO/commit/eeb3475252176002e26463dd93fb7d65e4cf00e7)) - wiltonloch
- moved yac calendar and datetime definition from CLEO to the reader python script - ([22efad5](https://github.com/yoctoyotta1024/CLEO/commit/22efad5edbd298355090ad4b8364a51023c2e9d8)) - wiltonloch
- update method to set library path - ([84c0669](https://github.com/yoctoyotta1024/CLEO/commit/84c06690a55e7d804d6c5c4c79875e0a8e62d5a0)) - clara.bayley

- - -

## [v0.20.0](https://github.com/yoctoyotta1024/CLEO/compare/8be2d7068d098a9a80c8abfe5b3aa64b86954ee7..v0.20.0) - 2024-06-21
#### Features
- add binary dir to cmake for CLEO externs - ([5cb5b52](https://github.com/yoctoyotta1024/CLEO/commit/5cb5b52c984d03ee6bbddfdd933473b2858df388)) - clara.bayley
#### Miscellaneous Chores
- format .gitignore - ([8be2d70](https://github.com/yoctoyotta1024/CLEO/commit/8be2d7068d098a9a80c8abfe5b3aa64b86954ee7)) - clara.bayley
#### Refactoring
- move src into roughpaper - ([102ad41](https://github.com/yoctoyotta1024/CLEO/commit/102ad41b8df53f16a2978ae6f51ed080642eb112)) - clara.bayley
- move src into roughpaper - ([1d4e3a0](https://github.com/yoctoyotta1024/CLEO/commit/1d4e3a02647ca101487979fc0c77a1cd2020b5ae)) - clara.bayley
- change compiler flags for correctness and debugging - ([0b5931b](https://github.com/yoctoyotta1024/CLEO/commit/0b5931b2f686f66402fcca3f5847cb72d8c5b709)) - clara.bayley
- move roughpaper into roughpaper/scratch - ([c657fd2](https://github.com/yoctoyotta1024/CLEO/commit/c657fd29ce22f0c51ab0fad8ed9ac4c4ddbf0dfd)) - clara.bayley
- move roughpaper into roughpaper/scratch - ([723f7d4](https://github.com/yoctoyotta1024/CLEO/commit/723f7d46991ab4d5605daeb29ccec3ae6ccb6e1b)) - clara.bayley

- - -

## [v0.19.0](https://github.com/yoctoyotta1024/CLEO/compare/eece61da0c37e039b8327af87d03cf558833297a..v0.19.0) - 2024-06-21
#### Features
- rough paper hpp to print size of gridboxes and superdroplets - ([94615ff](https://github.com/yoctoyotta1024/CLEO/commit/94615ffc9b5f98f8ee8aa5b28cdd1263c4b1d1de)) - clara.bayley
- rough paper hpp to print size of gridboxes and superdroplets - ([eece61d](https://github.com/yoctoyotta1024/CLEO/commit/eece61da0c37e039b8327af87d03cf558833297a)) - clara.bayley

- - -

## [v0.18.1](https://github.com/yoctoyotta1024/CLEO/compare/38b44b63d0bc5ba28187628564679bd8506b4063..v0.18.1) - 2024-06-21
#### Bug Fixes
- add fallthrough statements to show logic of switch is intentional - ([38b44b6](https://github.com/yoctoyotta1024/CLEO/commit/38b44b63d0bc5ba28187628564679bd8506b4063)) - clara.bayley
#### Refactoring
- more clarity and supress unused parameter warning in test build - ([27b6e25](https://github.com/yoctoyotta1024/CLEO/commit/27b6e25450f3bf7f97d286a58a70d102b4060f39)) - clara.bayley

- - -

## [v0.18.0](https://github.com/yoctoyotta1024/CLEO/compare/edb4dd3654440336ad52ab07f309532509014326..v0.18.0) - 2024-06-21
#### Features
- new helper script to run examples - ([5eb6ecb](https://github.com/yoctoyotta1024/CLEO/commit/5eb6ecb19ea10831ee72f71a31958470dfa0f9d7)) - clara.bayley
#### Refactoring
- move is_null_supers function content into collisions and delete redundant file - ([553cab0](https://github.com/yoctoyotta1024/CLEO/commit/553cab03ba4db8df44ba0bea785b2cf2cb23721a)) - clara.bayley
- move is_null_superdrop function to coalescence - ([ce99402](https://github.com/yoctoyotta1024/CLEO/commit/ce994023e9eb5e4ea257e39c3facf62ef50eb4de)) - clara.bayley
- minor doc fix and more explicit variable naming - ([1c89e09](https://github.com/yoctoyotta1024/CLEO/commit/1c89e0948d9bd4e39fcd928471c9413b0414bdd7)) - clara.bayley
- reduce minsubstep in speed test example - ([edb4dd3](https://github.com/yoctoyotta1024/CLEO/commit/edb4dd3654440336ad52ab07f309532509014326)) - clara.bayley

- - -

## [v0.17.0](https://github.com/yoctoyotta1024/CLEO/compare/4e9c3712c9acb976c71612cfe3fc59d40aae6d25..v0.17.0) - 2024-06-18
#### Bug Fixes
- fixed bug in substepping routine where incorrect rprev was used - ([9a9bf17](https://github.com/yoctoyotta1024/CLEO/commit/9a9bf174ab30a761fc937d9d0c107b9c21a5ee6c)) - clara.bayley
#### Continuous Integration
- added more restrictive compiler flags to treat all warnings as errors - ([9f6e01a](https://github.com/yoctoyotta1024/CLEO/commit/9f6e01ad633114969a5ceb5412348aa95de02a2a)) - wiltonloch
#### Documentation
- fix docstring typos - ([64b2daf](https://github.com/yoctoyotta1024/CLEO/commit/64b2daf0c0f59e69d3c3c79cd55edb1df36cd430)) - clara.bayley
#### Features
- new condensation method complete. removed unnecessary initial guess in substepping routine - ([f8034ef](https://github.com/yoctoyotta1024/CLEO/commit/f8034ef5b6f77708553fd1edb693df4365a2882d)) - clara.bayley
- adaptive sub-timestepping routine for solving non-unique g(Z) in implicit method - ([4cb6f2d](https://github.com/yoctoyotta1024/CLEO/commit/4cb6f2d4ae5d678b629f140868631f8e1a000ed1)) - clara.bayley
- adaptive sub-timestepping routine for solving non-unique g(Z) in implicit method - ([67fcf2b](https://github.com/yoctoyotta1024/CLEO/commit/67fcf2b1578b38111d264d18178ae75c52675069)) - clara.bayley
- add in tests of uniqueness to condensation ODE solving - ([153677e](https://github.com/yoctoyotta1024/CLEO/commit/153677e5f1f5541419ee3cfaa36c8f104676eb78)) - clara.bayley
#### Miscellaneous Chores
- **(yac_coupling)** added correct unsigned type to sizes and iterators - ([f5880ab](https://github.com/yoctoyotta1024/CLEO/commit/f5880abcf82ba094769d2af8972614d5c6fa6885)) - wiltonloch
- **(yac_coupling)** removing double declaration of coordinate bounds variables - ([63827e8](https://github.com/yoctoyotta1024/CLEO/commit/63827e8be0b126fc1bcbe3125911e0252735bfa0)) - wiltonloch
- formatting - ([f6d2c03](https://github.com/yoctoyotta1024/CLEO/commit/f6d2c03a6ad767602ba693c9ac151346883b184b)) - clara.bayley
- comment formatting - ([a94c0c4](https://github.com/yoctoyotta1024/CLEO/commit/a94c0c40d0aff901a5db4e273051c3c4180d0a81)) - clara.bayley
#### Refactoring
- nsupers correct for condensation testing setup - ([2df742f](https://github.com/yoctoyotta1024/CLEO/commit/2df742fef1e1f40cded9a5a8eec669f9973dd4c1)) - clara.bayley
- standard config for condensation - ([9f08cd4](https://github.com/yoctoyotta1024/CLEO/commit/9f08cd4093863bc15c4a6dc75a1e8c4901036e07)) - clara.bayley
- update config params for condensation's implicit method - ([d79a3c0](https://github.com/yoctoyotta1024/CLEO/commit/d79a3c07aef93aa9a14f3cdca562d252801e3a2f)) - clara.bayley
- reorganise functions for implicit euler method, work in progress - ([5a0e2e0](https://github.com/yoctoyotta1024/CLEO/commit/5a0e2e02e4c17b6184a293f864220ee92208f88d)) - clara.bayley
- using adia0D 0-D model setup for testing - ([cfdaa47](https://github.com/yoctoyotta1024/CLEO/commit/cfdaa47bbfde540e2e73a90adb5b8ad48927ede5)) - clara.bayley
- using adia0D 0-D model setup for testing - ([4e9c371](https://github.com/yoctoyotta1024/CLEO/commit/4e9c3712c9acb976c71612cfe3fc59d40aae6d25)) - clara.bayley

- - -

## [v0.16.1](https://github.com/yoctoyotta1024/CLEO/compare/a021794b0cd03474220ceccdc3d5eaab859a48d5..v0.16.1) - 2024-06-17
#### Bug Fixes
- fix missing re-scaling factor in phi calculation - ([a021794](https://github.com/yoctoyotta1024/CLEO/commit/a021794b0cd03474220ceccdc3d5eaab859a48d5)) - clara.bayley
#### Miscellaneous Chores
- docstring - ([89b52b4](https://github.com/yoctoyotta1024/CLEO/commit/89b52b455699bd6850a45c28c75e672b69dd3933)) - clara.bayley
- kokkos namespace explicit - ([a452daf](https://github.com/yoctoyotta1024/CLEO/commit/a452dafd2a4a8ee4bf8e2938f71773dfe25abb7f)) - clara.bayley

- - -

## [v0.16.0](https://github.com/yoctoyotta1024/CLEO/compare/aea9fd925dd71f85beede067e3a8992847d41ed4..v0.16.0) - 2024-06-17
#### Bug Fixes
- typos in executables - ([18f6072](https://github.com/yoctoyotta1024/CLEO/commit/18f60721d03160917e1091fc47a5933d5bbe38be)) - clara.bayley
- typo in build for CI - ([5ffc281](https://github.com/yoctoyotta1024/CLEO/commit/5ffc281532e1e4d5870e6fa1d7fd8219de9d8ded)) - clara.bayley
- set brekaup params in config construction (also minor formatting) - ([51559cd](https://github.com/yoctoyotta1024/CLEO/commit/51559cdc71bd1a838f9c11f478c43af44257c5ca)) - clara.bayley
#### Features
- add long kernel to breakup example - ([9b24c0c](https://github.com/yoctoyotta1024/CLEO/commit/9b24c0c9fa71ca19d6c1edd0b73e0d3ef814e2f5)) - clara.bayley
- more plots for breakup example - ([6471a14](https://github.com/yoctoyotta1024/CLEO/commit/6471a14e88597b4b1c2e480098336fa16677c5be)) - clara.bayley
- more plots for breakup example - ([9ac3f98](https://github.com/yoctoyotta1024/CLEO/commit/9ac3f98c36e62be411d1d83ebba2949b5ccb0473)) - clara.bayley
- include szakallurbich and testikstraub executables in breakup example - ([3999e55](https://github.com/yoctoyotta1024/CLEO/commit/3999e55948cae12bf0a28b8eebc629f3447f8dba)) - clara.bayley
- seperate cmake for boxmodelcollisions executable - ([b25719c](https://github.com/yoctoyotta1024/CLEO/commit/b25719cc8ff618d54f4bd0842aea24cc9bad529b)) - clara.bayley
- allow optional constnfrags nfrags a config parameter - ([0ad3329](https://github.com/yoctoyotta1024/CLEO/commit/0ad33298a0432f00ec05fbf608109c8c82286002)) - clara.bayley
- add new examples to CI build - ([3a48037](https://github.com/yoctoyotta1024/CLEO/commit/3a480377970fa76352afcbbc49001ecb0e38088a)) - clara.bayley
- include breakup in lowlist example - ([ab28c5e](https://github.com/yoctoyotta1024/CLEO/commit/ab28c5effb08dfcc5de3058d9df20c99c59bfec1)) - clara.bayley
- include breakup in lowlist example - ([c5a5dab](https://github.com/yoctoyotta1024/CLEO/commit/c5a5dab7d40db69eb4639deeced3129c11dd1dd4)) - clara.bayley
- new files for breakup executables - ([aea9fd9](https://github.com/yoctoyotta1024/CLEO/commit/aea9fd925dd71f85beede067e3a8992847d41ed4)) - clara.bayley
#### Miscellaneous Chores
- format - ([73f765f](https://github.com/yoctoyotta1024/CLEO/commit/73f765fcf13b0228e476e29bf6357761a7f763ba)) - clara.bayley
#### Performance Improvements
- beautify figures - ([54586f7](https://github.com/yoctoyotta1024/CLEO/commit/54586f7974d957b5ccfb1972a1f904356cf45a93)) - clara.bayley
#### Refactoring
- minor reordering of line styles - ([6b5113d](https://github.com/yoctoyotta1024/CLEO/commit/6b5113d1b58c1be15be23d251b9dc437b0b619e7)) - clara.bayley
- better plotting funcs - ([e8be503](https://github.com/yoctoyotta1024/CLEO/commit/e8be5030d599e6867d981529f74a1735974a5098)) - clara.bayley
- configure microphys for new breakup examples - ([6f9a436](https://github.com/yoctoyotta1024/CLEO/commit/6f9a4361c237075577eef57a48b873acf093409e)) - clara.bayley
- change nsupers in shima 2nd example - ([dc3fe62](https://github.com/yoctoyotta1024/CLEO/commit/dc3fe628b50632f06134868b51e3a1aae3873c74)) - clara.bayley

- - -

## [v0.15.0](https://github.com/yoctoyotta1024/CLEO/compare/f9fabde881fafaf459dbaaad846e72a3d38dd5f0..v0.15.0) - 2024-06-16
#### Documentation
- fixes to minor typos and examples docs - ([8422a27](https://github.com/yoctoyotta1024/CLEO/commit/8422a278c649b88828d5988d1449688433ad4abd)) - clara.bayley
#### Features
- new method for initial SD conditions for shima 2009 example - ([1dd7769](https://github.com/yoctoyotta1024/CLEO/commit/1dd7769f4893767ed516db0156d27a3bc69b5e0c)) - clara.bayley
#### Miscellaneous Chores
- fix formatting - ([9ab82ad](https://github.com/yoctoyotta1024/CLEO/commit/9ab82adbc02ec5feab79b3fa79b8c2009ab54522)) - clara.bayley
#### Performance Improvements
- configuration of shima 2009 example altered slightly - ([15fa197](https://github.com/yoctoyotta1024/CLEO/commit/15fa197594b6094ca50728cb47b9ea4e9dcd6d3b)) - clara.bayley
#### Refactoring
- split shima2009 examples from breakup example - ([f9fabde](https://github.com/yoctoyotta1024/CLEO/commit/f9fabde881fafaf459dbaaad846e72a3d38dd5f0)) - clara.bayley

- - -

## [v0.14.0](https://github.com/yoctoyotta1024/CLEO/compare/98566bb40ae938dbfd6065bdeb0efdfa46b68385..v0.14.0) - 2024-06-16
#### Features
- created subroutines to receive multiple horizontal levels from yac at once using yac collections - ([70cfa84](https://github.com/yoctoyotta1024/CLEO/commit/70cfa8406b964854b8df5f86f3525ee68b684e7f)) - wiltonloch
- changed cleo coupling to receive data from ICON's bubble test data reader - ([2bbf2bd](https://github.com/yoctoyotta1024/CLEO/commit/2bbf2bd9872ba789dc07211d51c4cff3a76251ad)) - wiltonloch
- added new python script to read data from ICON's bubble test - ([98566bb](https://github.com/yoctoyotta1024/CLEO/commit/98566bb40ae938dbfd6065bdeb0efdfa46b68385)) - wiltonloch
#### Miscellaneous Chores
- removed vertical slice-based subroutines to receive fields from yac - ([cc7ce27](https://github.com/yoctoyotta1024/CLEO/commit/cc7ce273b0200b0015358e218185a4219651f79f)) - wiltonloch
#### Performance Improvements
- added yac data containers as members of CartesianDynamics to have only one allocation and deallocation - ([6bed606](https://github.com/yoctoyotta1024/CLEO/commit/6bed6064a5f0146763fa5d7f0bc40bf1669724f0)) - wiltonloch
#### Refactoring
- redefined fields to receive all vertical levels at once - ([dbb572a](https://github.com/yoctoyotta1024/CLEO/commit/dbb572a555533f3a8c7e327f2e41ecda7df34a77)) - wiltonloch

- - -

## [v0.13.0](https://github.com/yoctoyotta1024/CLEO/compare/b8cdcab0aedab3e2f0f2e0fc9086a06f0fe566d6..v0.13.0) - 2024-06-11
#### Features
- optional lat lon bound parameters applied to vertex coordinates creation - ([cb88b57](https://github.com/yoctoyotta1024/CLEO/commit/cb88b576856d577c621e20aa4bb4d618efd0d50f)) - wiltonloch
- added lat and lon bounds as optional parameters to yac dynamics - ([47f5b77](https://github.com/yoctoyotta1024/CLEO/commit/47f5b77391a90f499e50108240b469be9cc29d88)) - wiltonloch
#### Miscellaneous Chores
- renamed python data reader specifying target cleo data - ([05fef44](https://github.com/yoctoyotta1024/CLEO/commit/05fef44fb967d8cb0984bd6630c240a7e4240b88)) - wiltonloch
#### Refactoring
- added explit lat and lon bounds for vertex coordinates calculations - ([869854a](https://github.com/yoctoyotta1024/CLEO/commit/869854a5aa6523dd2f2355a302f8712ba87a7e4e)) - wiltonloch
- encapsulated vertex coordinates calculation - ([751e426](https://github.com/yoctoyotta1024/CLEO/commit/751e426093b0bd329d0836888fd55ed1335fb651)) - wiltonloch
- made cell center lat and lon calculations dependent only on vertex lat and lon - ([5364e6a](https://github.com/yoctoyotta1024/CLEO/commit/5364e6a6e21d2467b06fe6dee6eced3a78332909)) - wiltonloch
- encapsulated grid and points definition and removed vertex lat and lon arrays from CartesianDynamics - ([4c9d488](https://github.com/yoctoyotta1024/CLEO/commit/4c9d488fadc8200019ca7449fbc5f0587d92eabe)) - wiltonloch
- grouped yac id declarations and replaced total cells and edges by yac calls for target arrays creation - ([b8cdcab](https://github.com/yoctoyotta1024/CLEO/commit/b8cdcab0aedab3e2f0f2e0fc9086a06f0fe566d6)) - wiltonloch

- - -

## [v0.12.0](https://github.com/yoctoyotta1024/CLEO/compare/4a7710eb552ff14fc2caa21baf712b2d86ecf663..v0.12.0) - 2024-06-06
#### Features
- new bash script to run yac divfree2d example - ([4a7710e](https://github.com/yoctoyotta1024/CLEO/commit/4a7710eb552ff14fc2caa21baf712b2d86ecf663)) - clara.bayley
#### Miscellaneous Chores
- add notes - ([84160fd](https://github.com/yoctoyotta1024/CLEO/commit/84160fdbc191108755d1d60904b664454e8ab41d)) - clara.bayley

- - -

## [v0.11.0](https://github.com/yoctoyotta1024/CLEO/compare/5ddd84535e29462079c6d61a60fcdd3c4bfcd3e3..v0.11.0) - 2024-06-06
#### Bug Fixes
- corrected args passed to run speedtest example - ([fa5e4c5](https://github.com/yoctoyotta1024/CLEO/commit/fa5e4c5ba6583eaeb2cc62a451528b31e0c47121)) - clara.bayley
- remove unwanted single thread commands - ([648045a](https://github.com/yoctoyotta1024/CLEO/commit/648045ab2add330f3770d589aa7d75cfcad085d5)) - clara.bayley
- missing pass by reference for xzarr - ([53a72fc](https://github.com/yoctoyotta1024/CLEO/commit/53a72fc0e4ab98ba51cf75681fa1001c89db3539)) - clara.bayley
- minor bug fixes - ([2e956f6](https://github.com/yoctoyotta1024/CLEO/commit/2e956f6b53773d7822826bde467cb79181b7817c)) - clara.bayley
- sdmmonitor concept cannot take d_gbxs in constraints - ([51cbb09](https://github.com/yoctoyotta1024/CLEO/commit/51cbb09a87f39cfc71800ef70dd186dbecefb861)) - clara.bayley
- minor bug fixes e.g. capture pointer to class instance in lambdas - ([30567a7](https://github.com/yoctoyotta1024/CLEO/commit/30567a7eebda686d81f457bea9b25fbd52434d7a)) - clara.bayley
- single thread barrier to monitor condensation - ([5fd86bb](https://github.com/yoctoyotta1024/CLEO/commit/5fd86bb33d530985f6a534efac1f8dd3d2e88df4)) - clara.bayley
- count is modifyable in view - ([b95b6f2](https://github.com/yoctoyotta1024/CLEO/commit/b95b6f25917a7ba661bb56154be5a39dba319388)) - clara.bayley
#### Documentation
- note on missing constraint on SDMMonitor - ([020658b](https://github.com/yoctoyotta1024/CLEO/commit/020658b3cf4f9d030f32515e4695c7dc618f4e27)) - clara.bayley
#### Features
- new files for monitoring mass moments during SDM - ([5919ca9](https://github.com/yoctoyotta1024/CLEO/commit/5919ca9ac8dcc429978f16eaea9e33c7adacb3cc)) - clara.bayley
- new files for monitoring mass moments during SDM - ([5ddd845](https://github.com/yoctoyotta1024/CLEO/commit/5ddd84535e29462079c6d61a60fcdd3c4bfcd3e3)) - clara.bayley
#### Miscellaneous Chores
- renaming and remove comments - ([e5ea219](https://github.com/yoctoyotta1024/CLEO/commit/e5ea2191ed23662a3222330ac3ec1ec0a28f28c7)) - clara.bayley
#### Refactoring
- extend time for SLURM jobs - ([a8b9494](https://github.com/yoctoyotta1024/CLEO/commit/a8b94948ce858a1331996049ad05a49376f0d77b)) - clara.bayley
- flag in bash to check valid buildtype - ([57e137b](https://github.com/yoctoyotta1024/CLEO/commit/57e137b53eea998587b87826bddbf5491ce96944)) - clara.bayley
- delete attempt to set LAPACK flags in sundials - ([c61e489](https://github.com/yoctoyotta1024/CLEO/commit/c61e4891158698d00f2298f7419cdc3401d92cf1)) - clara.bayley
- rename and use actual calc massmoms function - ([511f767](https://github.com/yoctoyotta1024/CLEO/commit/511f7672dd3adf7f5f00a9b9716dfd4a13d90d85)) - clara.bayley
- move rainmassmoms calc out of functor into seperate function - ([53a7f39](https://github.com/yoctoyotta1024/CLEO/commit/53a7f393b5ac06f438f8150b62feb29ee4f49c28)) - clara.bayley
- include func concept in type - ([460c1c6](https://github.com/yoctoyotta1024/CLEO/commit/460c1c66e3825e1a28e3bbcb90e4621274338dcb)) - clara.bayley
- include condensation monitor in eurec4a example - ([59b61b2](https://github.com/yoctoyotta1024/CLEO/commit/59b61b2d9c2a8e56016e349ddfa3fa33cfbe03a9)) - clara.bayley
- single thread for calc massmoments call in monitoring motion - ([2087fde](https://github.com/yoctoyotta1024/CLEO/commit/2087fdeebb77265a7d40e813c90c509040cb01e2)) - clara.bayley
- place plugs for monitoring condensation, microphysics and motion - ([d3e7f02](https://github.com/yoctoyotta1024/CLEO/commit/d3e7f0238b82657caab1969ea71bf899b818c787)) - clara.bayley
- place plugs for monitoring condensation, microphysics and motion - ([09e24a0](https://github.com/yoctoyotta1024/CLEO/commit/09e24a0e429290314afd4c46e1e41f6dca31b7e3)) - clara.bayley
- move function for mass moments calc out of functor - ([cff8c3a](https://github.com/yoctoyotta1024/CLEO/commit/cff8c3af818d00a4a28e6bc4a6eb3abaa2ca3188)) - clara.bayley
- extend mass moments monitors to all 3 mass moments for microphys and motion - ([865c8cf](https://github.com/yoctoyotta1024/CLEO/commit/865c8cf130a8e21c57dbe0753534842850158f93)) - clara.bayley
- move observers array creation to seperate file - ([6dc11a9](https://github.com/yoctoyotta1024/CLEO/commit/6dc11a9265d7e17ce071deaa771e77ad2a9fa76b)) - clara.bayley
- rename and move mass moments observer to seperate file - ([e532a0c](https://github.com/yoctoyotta1024/CLEO/commit/e532a0c58fa8cd288d6f87e8c28ad5c9e2ea0b8e)) - clara.bayley
- rename and move mass moments observer to seperate file - ([5486919](https://github.com/yoctoyotta1024/CLEO/commit/5486919a0b2e9e254c8a077b27236b4825014b40)) - clara.bayley
- rename and move mass moments observer to seperate file - ([9e6f8bf](https://github.com/yoctoyotta1024/CLEO/commit/9e6f8bf566f642df1201c1b134708bf81d1c53bc)) - clara.bayley
- include mass moments monitoring observer - ([ee74de9](https://github.com/yoctoyotta1024/CLEO/commit/ee74de9b5e1ff02aae9e3f876692a3369f42341f)) - clara.bayley
- extend SDMMonitor concept for mass moments monitoring - ([f172672](https://github.com/yoctoyotta1024/CLEO/commit/f17267215c0440330c32f754344bf288cb22ec27)) - clara.bayley

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
