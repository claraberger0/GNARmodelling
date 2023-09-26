
# GNAR modelling for congestion costs in ERCOT

<h3 align="center">Using network autoregressive models to predict nodal congestion costs in the ERCOT electrical grid</h3>
<!-- PROJECT LOGO -->
<br />
<div align="center">
  <a href="https://github.com/claraberger0/GNARmodelling">
    <img src="figures/tx_network_.png" alt="Logo" width="500" height="500">
  </a>

</div>

This project was carried out in collaboration with [[Aurora Energy Research](https://auroraer.com/)].

<!-- ABOUT THE PROJECT -->
## About The Project

Codebase for Clara Berger's MSc Dissertation at University of Oxford: Using network autoregressive models to predict nodal congestion costs in the ERCOT electrical grid

Abstract: 
Nodal pricing models require that energy prices be modelled at a number of in- teracting nodes simultaneously. In this report, we model the congestion cost com- ponent of energy prices at various nodes in the Electric Reliability Council of Texas (ERCOT) grid. We utilise Generalised Network Autoregressive (GNAR) models to constrain the interactions between various nodes based on network structure. Our method contributes a novel approach in electricity price modelling by employing an autoregressive model at the nodal level that incorporates grid structure. The models perform well in fitting the nodal historical congestion cost data, but display limited forecasting capabilities when predicting more than one hour in the future. Assum- ing the edges are present a priori leads to a robust method of dimension reduction for Vector Autoregressive modelling, but limits us to only modelling between 8 and 59 nodes at a time. Overall, incorporating neighbour effects improves the models’ prediction capabilities over a univariate AR model for 11 out of 16 combinations of time series and network structure.

<!-- GETTING STARTED -->
## Getting Started

Key packages required:
R:
* `GNAR`  
* `igraph`
* `xts`

Python:
* `geopandas`
* `pandas`

<p align="right">(<a href="#readme-top">back to top</a>)</p>



<!-- ROADMAP -->
## Contents

### code
  - `run_get_grid_geometry.py` python file for constructing the electric grid of Texas using public data, find all node locations to create an adjacency matrix as well merge nodes from the public data with those from Aurora
  - `get_grid_geometry.py` python file containing the functions for `run_get_grid_geometry.py`

  - `GNAR.R` R file that creates a network, manipulates the historical congestion cost data, fits the GNAR model, assess model performance, and forecasts future congestion costs based on the fit model
  - `run_GNAR.R` R file containing the functions for `run_get_grid_geometry.py`
    
  - `network_statistics.R` R file that creates an `igraph` network for the Texas electric grid and then extracts key network statistics

### network_comp
  - `network_comparison.R` R file that compares the networks created using public data and the actual electric grid of Texas

### grid_shape
  - `us_grid.shp` shapefile from the feature layer U.S. Electric Power Transmission Lines based on data from Homeland Infrastructure Foundation-Level Data (HIFLD) [
U.S. Electric Power Transmission Lines]([https://www.census.gov/geographies/mapping-files/time-series/](https://www.arcgis.com/home/item.html?id=d4090758322c4d32a4cd002ffaa0aa12&view=list&sortOrder=desc&sortField=defaultFSOrder))

### cb_2022_us_state_5m
  - `cb_2022_us_state_5m.shp` shapefile of US states from [US Census TIGER API](https://www.census.gov/geographies/mapping-files/time-series/) used to cut out the transmission lines in Texas

### figures
  - `tx_network_.png` figure of the network created from the Texas elextric grid with transmission lines as edges based on public data


<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- CONTACT -->
## Contact

Clara Berger - cberger4@wellesley.edu

Project Link: [https://github.com/claraberger0/GNARmodelling](https://github.com/claraberger0/GNARmodelling)

<p align="right">(<a href="#readme-top">back to top</a>)</p>





<!-- MARKDOWN LINKS & IMAGES -->
<!-- https://www.markdownguide.org/basic-syntax/#reference-style-links -->
[contributors-shield]: https://img.shields.io/github/contributors/claraberger0/GNARmodelling.svg?style=for-the-badge
[contributors-url]: https://github.com/claraberger0/GNARmodelling/graphs/contributors
[forks-shield]: https://img.shields.io/github/forks/claraberger0/GNARmodelling.svg?style=for-the-badge
[forks-url]: https://github.com/claraberger0/GNARmodelling/network/members
[stars-shield]: https://img.shields.io/github/stars/claraberger0/GNARmodelling.svg?style=for-the-badge
[stars-url]: https://github.com/claraberger0/GNARmodelling/stargazers
[issues-shield]: https://img.shields.io/github/issues/claraberger0/GNARmodelling.svg?style=for-the-badge
[issues-url]: https://github.com/claraberger0/GNARmodelling/issues
[license-shield]: https://img.shields.io/github/license/claraberger0/GNARmodelling.svg?style=for-the-badge
[license-url]: https://github.com/claraberger0/GNARmodelling/blob/master/LICENSE.txt
[linkedin-shield]: https://img.shields.io/badge/-LinkedIn-black.svg?style=for-the-badge&logo=linkedin&colorB=555
[linkedin-url]: https://linkedin.com/in/linkedin_username
[product-screenshot]: images/screenshot.png
[Next.js]: https://img.shields.io/badge/next.js-000000?style=for-the-badge&logo=nextdotjs&logoColor=white
[Next-url]: https://nextjs.org/
[React.js]: https://img.shields.io/badge/React-20232A?style=for-the-badge&logo=react&logoColor=61DAFB
[React-url]: https://reactjs.org/
[Vue.js]: https://img.shields.io/badge/Vue.js-35495E?style=for-the-badge&logo=vuedotjs&logoColor=4FC08D
[Vue-url]: https://vuejs.org/
[Angular.io]: https://img.shields.io/badge/Angular-DD0031?style=for-the-badge&logo=angular&logoColor=white
[Angular-url]: https://angular.io/
[Svelte.dev]: https://img.shields.io/badge/Svelte-4A4A55?style=for-the-badge&logo=svelte&logoColor=FF3E00
[Svelte-url]: https://svelte.dev/
[Laravel.com]: https://img.shields.io/badge/Laravel-FF2D20?style=for-the-badge&logo=laravel&logoColor=white
[Laravel-url]: https://laravel.com
[Bootstrap.com]: https://img.shields.io/badge/Bootstrap-563D7C?style=for-the-badge&logo=bootstrap&logoColor=white
[Bootstrap-url]: https://getbootstrap.com
[JQuery.com]: https://img.shields.io/badge/jQuery-0769AD?style=for-the-badge&logo=jquery&logoColor=white
[JQuery-url]: https://jquery.com 

