# Changelog

## v1.50.4

* Feature: Added a way to convert a `fotf` object to a symbolic expression. The new function is `fotf2expr`. Note that this feature requires the **Symbolic Math Toolbox**. There are some updates to existing functions as well. For example, the `poly2str` function now supports a 4th argument which is the kind of multiplication symbol to be inserted according to need -- by default, it is empty, i.e. no multiplication symbol is inserted, and this is used in FOMCON GUIs when you see the string expressions for FOTF numerator and denominator.
- Bugfix: **[breaking]** FOMCON function `poly2str` being shadowed by a Control System toolbox function of the same name. The solution is to keep this function in FOMCON for backwards-compatibility as an alias to `fpoly2str`. To keep consistency, the function `str2poly` has also been recreated as `str2fpoly`.

- Bugfix: (regression?) simulating FOTF models with a delay term resulted in an incorrect concatenation of sample vectors.

## v1.50.3

* Resolved regression bug with the release whereby all the necessary files and folders were not added to MATLAB path upon installing the toolbox.

* Slightly updated the documentation in the README file concerning the new installation method.

* Added the project file for packaging the toolbox.

## v1.50.2

Hotfix for MATLAB R2022a wherein an error is produced when MATLAB is started due to a missing line in the bundled `info.xml` file. The file has been removed for now, but may return in the future, when the documentation is updated.

## v1.50.1

Streamlined releases through the Mathworks File Exchange/GitHub pipeline. There were no changes to the toolbox code.

## v1.50.0

Major changes:

* Initial publication of the changelog. 
* Update to the Grünwald-Letnikov based FODE solver: added higher-order numerical solver options -- the default and only available option up to this point was the first order solver; now, the user gets the option of also choosing second and third order numerical schemes for computing the time domain response of systems described by fractional-order transfer function objects. The setting can be configured in the new option `GL_Order`  found under `Simulations` in the toolbox configuration.
* Completely rewritten the configuration options handler. This is, in some ways, **a breaking change**, but at the same time one that should not affect the majority of users. Whereas before the user had to store or recall toolbox-wide configuration options him or herself and the configuration was stored in MATLAB workspace rather than in the user's file system, the new approach is to: (1) Store the configuration as a file in the user's home directory under the `.fomcon-matlab` folder. This way, the configuration options are persistent from session to session. (2) Perform a check on the options, so if new options are added in the future, there is no backwards compatibility problem. On the other hand, the validity of the inserted values/ranges for the configuration parameters are currently not checked, this is a feature that was also missing previously and will be added in the future.
* Added a new demo which pertains to the above mentioned solver update to show the differences in the solutions that can be obtained using various methods and with different solver orders.
* Released a packaged toolbox version of the toolbox.

Minor changes:

* Some minor bugs/compatibility problems have been fixed.

Known issues:

* In newer versions of MATLAB, a warning is displayed when the UI of the toolbox configuration dialog is displayed about the deprecation of certain Java components. In FOMCON, a contribution called propertiesGUI is used to display this configuration box and it has not been updated since 2015. In a future version of the toolbox, it is likely to be replaced with a custom, programmatic UI.

## v1.22.0

This can be considered the fundamental version of the toolbox published in 2015. There is no changelog prior to or for this version, but the development of the toolbox up to 2015 can be traced by following the publications of A. Tepljakov et al. from 2011 starting with the first [journal publication](https://ijmcs.dmcs.pl/documents/10630/20046/JMCS_2_2011-3.pdf) in which the initial version of the toolbox was introduced.
