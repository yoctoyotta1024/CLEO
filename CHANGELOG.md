# Changelog
All notable changes to this project will be documented in this file. See [conventional commits](https://www.conventionalcommits.org/) for commit guidelines.

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
