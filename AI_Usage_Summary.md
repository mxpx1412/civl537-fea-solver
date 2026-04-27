# AI Usage Summary

## Overview

### AI Tools Used

Claude was the only AI consulted throughout this project. 

## What was Written by Me vs. AI-Assisted

To maximize personal learning, AI use was limited to specific cases. A list of
AI assisted items are listed
[below](#list-of-issues-resolved-with-ai-assistance). In general, the following
content was produced by myself (excluding content already provided by course):

### Contents by Myself

+ `README.md`: typed own write up (referencing sources and providing
  bibliography for referenced contents), created own figures using `draw.io`.
  Ai checked select parts (see below) but no content is copied from AI.
+ `src/*` modules: implemented own solutions, AI checked or provided advise on
  points of confusion where applicable (see below). So original work include
  (aside from already provided code):
    + `elements.py`
    + `mesh.py`
    + `assembly.py`
    + `solver.py`
    + `postprocess.py`
+ Other files are largely already provided, but changes are made as appropriate
  as outlined in `README`

### List of Issues Resolved with AI Assistance

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
+ **Used Claude to help with final app deployment**, uses:
    + Asked Claude for advice when `app` updates were not reflected in webapp.
    + Asked for suggestion on UI tweaks.
    + Asked for advice on filtering for "near hole" elements, since initial
      code was not working. Implemented selecting elements by node approach.
    + Asked for advice on fixing deformation plots disappearing when tweaking
      deformation scale.
+ **Used Claude to better understand locking**
    + Asked Claude for high level overview.
    + Asked Claude to provide sources for further reading. 
    + Read and provided explanation in own words in README.

## One Case Where AI Was Wrong

Due to intentionally limiting AI usage, there were limited AI errors detected.
However, one prominent example occurred during the `app.py` deployment
(transition from starter to final). I found that my changes were not being
reflected in the app (both local and webapp). I asked Claude for advice,
however, Claude provided an incorrect answer, where it initially insisted my
shell/terminal is the problem. Claude claimed my docker container was stuck due
to shell variables set by my terminal preventing an update. However, upon
further investigation, I found that the docker image is referencing my old
`app_starter.py` file because it cached the old file, so it was not referencing
the updated file and that was preventing change. I asserted this finding to
Claude, and it finally gave the correct advise of rebuilding the docker image
and clearing the cache. 

## What I Learned

Overall, AI like Claude is very useful in assisting targeted problems. However,
I do not think I would learn nearly as much about the FEA process if I asked it
to create everything for me. Further, the mistakes that AI makes can be "silent
but deadly". Without independent judgment, know-how and the ability to assert
one's own expertise, we can very easily be misled by AI to go on wild
"goose-chases" that do not solve the problem. Moving forward, I will continue
to use AI assistant in targeted problems, but ensure that I am gaining my own
understanding and reserve the ability to exercise independent judgemnt, as to
prevent AI mistakes from slipping into important work.
