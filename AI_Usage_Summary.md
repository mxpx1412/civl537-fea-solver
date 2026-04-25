# AI Usage Summary

## Overview

AI assistance was employed in the course of this project as outlined below.

## List of Issues Resolved with AI Assistance

+ **Used Claude to assist with initial deployment**:
    + *Issue*: Streamlit deployment from repo initially failed.
    + *Fix*: 
        + Identified issue: mismatched Streamlit Python version (Streamlit
          defaults to 3.14, 3.11 required).
        + Fixed Python version setting in Streamlit, deployment succeeded.
+ **Used Claude to help with README formatting**: 
    + *Issue*: Markdown has limited image positioning capabilities, images looked off. 
    + *Fix*: adopted Claude suggestion to use HTML to position pictures.
+ **Used Claude to assist with error message from pytest**:
    + *Issue*: Failed to run test
        + Ran: `docker compose run fea-solver pytest tests/test_elements.py -v`
        + Error: `E   ModuleNotFoundError: No module named 'src'`
    + *Fix*: adopted Claude fix by adding `pyproject.toml` and adding root path
      for `pytest` initiation
+ **Used Claude to check load vector algebra**:
    + *Issue*: `R.sum` tests were failing (does not sum to `P`)
    + *Fix*: uploaded `README` to Claude and asked for algebra check,
      identified sign issue. Fixed algebra sign and updated `assembly.py` to
      match, test passed after
+ **Used Claude for miscellaneous Q&A on Python / Python libraries**, examples:
    + Asked Claude to explain `lil_matrix` and `csr` matrix implmentation
    + Asked Claude the best way to select certain rows/columns of matrices,
      received recommendation using `np.ix_`, received
+ **Used Claude to identify failed `test_solver` test** reason:
    + *Issue*: cantilever deflection limit failed in `test_solver` test
    + *Fix*:
        + First independently checked, suspected issue with computation of `I`
          and use of `0.01` in thickness in the test.
        + Uploaded `src` and `test` files to Claude and asked Claude to confirm
          if this is the issue or if there are other bugs.
        + Claude confirmed that the thickness is the main issue. Modified test
          and added docstring after Claude confirmation.
        + Claude also identified side issue of indexing in
          `assemble_R_uniform_tension`, fixed per Claude suggestion.

