pro print_moments,vector,mom,silent=silent

  if n_params() eq 0 then begin
      print,'-syntax print_moments,vector,silent=silent'
      print,'print out moments of vector and store them in mom'
      return
  endif
  
  mom=moment(vector)
  
  if (not keyword_set(silent) ) then begin
      print,'Mean: ',mom[0]
      print,'Variance: ',mom[1]
      print,'Standard Deviation: ', sqrt(mom[1])
      print,'Skewness: ',mom[2]
      print,'Kurtosis: ',mom[3]
  endif
  
  return
end
