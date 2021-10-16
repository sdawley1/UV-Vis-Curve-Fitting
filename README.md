# UV-Vis-Curve-Fitting
AS.030.421 Homework #5: UV/Vis Curve Fitting

This was the definitely the most challenging assignment yet and I do not at all expect them to get easier. I'm also not convinced that I arrived at the correct result, but I had fun coding this nonetheless. If I was (at least close to being) correct, this is probably not the method that was used for the actual research project.

  For context, here is the publication that laid the groundwork for why we actually wrote this program in the first place: https://pubs.acs.org/doi/full/10.1021/la0351764

  The point of the research was to elucidate a method of determining the sample composition of a surfactant using UV/Vis data (there's more to it than that - read the paper if you're interested). I would certianly never have thought anyone could determine sample composition from UV spectroscopy but then again this work was done by people a lot smarter than me. Ultimately, this method is not only incredibly interesting but also it offers some interesting insight into the interrelatdness of mathematics and chemistry. 
   
  The program I wrote relies on fitting UV/Vis data using a gaussian distribution and creating a linear program from a system of distributions which relies on the physical restrictions in place. The coding aspect of this algorithm was fairly straightforward, all things considered. The mathematical foundation behind why any of it works is also pretty cool and relies on some nice properties of the gaussian distribution. I mention them in the attached Jupyter notebook. 

UPDATE: This approach to the problem did in fact arrive at the correct result. In fact, this was the same method used when the research was originally conducted.
