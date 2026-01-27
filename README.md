# [Ontogenetic shifts in wood anatomy and leaf traits in tropical dry forests](https://nph.onlinelibrary.wiley.com/doi/10.1111/nph.70725)

### New Phytologist

### Peter J. Williams, Elise F. Zipkin, Andrés González-Melo, Beatriz Salgado-Negret, Roy González-M., Natalia Norden, Juan Pablo Benavides-Tocarruncho, Juan Manuel Cely, Julio Abad Ferrer, Daniel García-Villalobos, Fabián Garzón, Álvaro Idárraga-Piedrahita, René López-Camacho, Esteban Moreno, Jhon Nieto, Camila Pizano, Juliana Puentes-Marín, Nancy Pulido, Katherine Rivera, Felipe Rojas-Bautista, Viviana Salinas, Juan Felipe Solorzano, María Natalia Umaña

### Data DOI: [https://doi.org/10.7302/kjr3-d129](https://doi.org/10.7302/kjr3-d129)

### Please contact the first author for questions about the code or data: Peter J. Williams (peter.j.williams.110@gmail.com or will4047@msu.edu)

---------------------------------

## Abstract
- Trees experience changes in environmental conditions and biophysical constraints as they grow that may lead to shifts in functional traits. While ontogenetic shifts in leaf traits are relatively well understood, changes in wood anatomy traits from seedlings to adults are less clear, especially in tropical dry forests where drought strongly influences adult wood traits.
- We collected wood and leaf functional trait data for seedlings and adult trees in four dry forest sites along a rainfall gradient in Colombia. Using a Bayesian framework, we quantified intraspecific trait shifts between seedlings and adults, and we compared functional axes and trait correlations between these stages.
- Leaf traits shifted along a clear functional axis from more acquisitive seedlings to more conservative adults, but changes in vessel traits were more complicated. Vessel diameter increased across ontogeny as plant height increased. Wood anatomy traits were less tightly coupled for seedlings than for adults, such that the hydraulic safety–efficiency trade-off observed in adults appeared to be weaker.
- Our results indicate that wood anatomy traits do not show coordinated ontogenetic shifts among traits because traits associated with hydraulic safety and efficiency are less integrated in seedlings. In dry forests, hydraulic trade-offs become stronger as trees grow taller. 

## [Code](Code)
- **[1_TDF_ontogeny_trait_shifts.R](Code/1_TDF_ontogeny_trait_shifts.R)**: Analyze ontogenetic trait shifts from seedlings to adults, analyze traits shifts across sites (Results: Questions 1 & 2).
- **[1_TDF_ontogeny_PCA.R](Code/1_TDF_ontogeny_PCA.R)**: Calculate Bayesian PCA with varimax rotation to calculate functional trait axes (Results: Question 3).
- **[1_TDF_ontogeny_trait_correlations.R](Code/1_TDF_ontogeny_trait_correlations.R)**: Calculate trait correlations for adults and for seedlings, calculate differences in slopes and intercepts between adults and seedlings (Results: Question 3).
- **[2_TDF_ontogeny_figures.R](Code/1_TDF_ontogeny_figures.R)**: Create all figures (Figs. 1-4 & Fig. S1).

## Data
- The file "traits_data.csv" is available in the Deep Blue Data repository at: https://doi.org/10.7302/kjr3-d129.
