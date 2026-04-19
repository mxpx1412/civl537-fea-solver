# AI Usage Summary

AI assistance was employed in the course of this project as outlined below.

## Overview


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
