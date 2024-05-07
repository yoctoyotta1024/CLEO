# Changelog
All notable changes to this project will be documented in this file. See [conventional commits](https://www.conventionalcommits.org/) for commit guidelines.

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
