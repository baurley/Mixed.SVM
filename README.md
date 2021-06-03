![alt text](http://i.imgur.com/NB0Y8Su.png "BioRealm Logo")

# Mixed.SVM

A Bayesian mixed effects Support Vector Machine for learning and predicting daily alcohol use disorder patterns.

## Table of Contents

- [**Project Organization**](#project-organization)
- [**Team Members**](#team-members)
- [**Project Motivation**](#project-motivation)
- [**How Do I Run It/How Does It Work**](#how-do-i-run-it-how-does-it-work)
- [**Roadmap**](#roadmap)
- [**Installation**](#installation)
- [**Contributing**](#contributing)
- [**Contributors**](#contributors)
- [**Changelog**](#changelog)
- [**Security**](#security)
- [**License**](#license)

# <a name="project-organization"></a>Project Organization

    ├── AUTHORS
    ├── CHANGELOG.md
    ├── CODE_OF_CONDUCT.md
    ├── CONTRIBUTING.md
    ├── CONTRIBUTORS
    ├── data
    │   └── processed
    │       └── example.RData                    <- simulated daily substance use data
    ├── docs
    ├── LICENSE
    ├── notebooks                                <- Jupyter notebooks
    │   ├── daily-drinking-notebook-demo.ipynb  <- interactive analysis of simulated daily substance use data
    │   └── daily-drinking-notebook.ipynb        <- analysis of the ABQ DrinQ data
    ├── README.md
    ├── reports                      
    │   ├── figures
    │   └── results.RData                    <- saved results from running mixed-svm
    ├── SECURITY
    ├── src                            
    │   └── R
    │       └── mixed-svm.R                    <- source code for the MCMC algorithm and summaries
    └── VERSION

---

# <a name="team-members"></a>Team Members

- "James W. Baurley" <baurley@biorealm.ai>
- "Christopher S. McMahan" mcmaha2@clemson.edu

# <a name="project-motivation"></a>Project Motivation

Alcohol use disorder (AUD) is a heterogeneous disorder. Classification algorithms may accelerate AUD research by parsing heterogeneity in clinically meaningful ways. Machine learning algorithms are capable of analyzing high dimensional data, but lack the ability to model the more complicated features inherent in AUD data. Inspired by a study of heavy drinkers that collected longitudinal neuroimaging, daily logs, and extensive self-report and interview assessments (ABQ DrinQ), we developed an algorithm for learning and predicting AUD patterns. We recast support vector machines (SVMs), a common classification technique, into a Bayesian model extended to handle mixed effects. We then apply these methods to ABQ DrinQ to develop models of alcohol use patterns and risk factors for heavy drinking. We identified gender, age, and time varying usage of nicotine, cannabis, and other drugs as significant risk factors of heavy drinking behavior. The model classifies 84% of heavy drinking days correctly. The algorithms and tools to summarize the results are packaged into a source code repository for researchers to explore, along with examples in an interactive Jupyter notebook. AUD researchers may use this toolkit to characterize daily use patterns and identity risk factors of daily use in their own studies. Ultimately, understanding patterns of alcohol use and risk for alcohol use could be used for developing individualized interventions. 

# <a name="how-do-i-run-it-how-does-it-work"></a>How Do I Run It/How Does It Work

Open the `daily-drinking-notebook.ipynb` notebook in Jupyter to follow the analysis presented in the paper. Open `daily-drinking-notebook-demo.ipynb` to interactively execute the analysis on a provided simulated dataset of daily substance use.



Please see the paper for details on the methodology.

# <a name="roadmap"></a>Roadmap

We plan to continue to develop the toolkit to help alcohol use disorder (AUD) researchers gain new insights from their study data.  The toolkit consists of machine learning algorithms adapted to handle mixed effects (e.g., Support Vector Machine) and Jupyter notebooks to help the user understand and apply the algorithms.

# <a name="installation"></a>Installation

Please follow the Anaconda documentation: [Using the R programming language in Jupyter Notebook &#8212; Anaconda documentation](https://docs.anaconda.com/anaconda/navigator/tutorials/r-lang/)



If the fonts for figures are not displaying correctly, please install the `r-cairo` package in Anaconda Navigator. 

# <a name="contributing"></a>Contributing

 Please see the <a href="CONTRIBUTING.md">`CONTRIBUTING.md`</a> file in the project's root directory for detailed information on how to contribute to this project, including detailed code, documentation, and Git style guidelines.

# <a name="contributors"></a>Contributors

For a complete list of contributors, please see the <a href="CONTRIBUTORS">`CONTRIBUTORS`</a> file in the project's root directory.

# <a name="changelog"></a>Changelog

Please see the <a href="CHANGELOG">`CHANGELOG`</a> file in the project's root directory.

# <a name="security"></a>Security

Please see the <a href="SECURITY">`SECURITY`</a> file in the project's root directory for information on reporting any suspected security vulnerabilities.

# <a name="license"></a>License

Please see the <a href="LICENSE">`LICENSE`</a> file in the project's root directory.
