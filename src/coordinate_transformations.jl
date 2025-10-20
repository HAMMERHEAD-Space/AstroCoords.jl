# I cannot see src/coordinate_transformations.jl

# I need to see the existing file structure to properly extend the COE transformation functions
# without breaking existing elliptical orbit code. Please provide the current contents of:
# - src/coordinate_transformations.jl
# - Any related files that define Vec3 types or helper functions
# - The module structure in src/AstroCoords.jl (to understand exports)

# Once I can see the existing implementation, I will:
# 1. Extend coe_to_rv() to handle e≥1 cases
# 2. Extend rv_to_coe() to handle e≥1 cases  
# 3. Add mean_to_true_anomaly() and true_to_mean_anomaly() helpers
# 4. Preserve all existing functionality for elliptical orbits
# 5. Add comprehensive docstrings citing Vallado 2013 equations