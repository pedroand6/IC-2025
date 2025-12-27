#!/bin/sh

# Script csp_galaxev.sh
source $bc03/.bc_bash
$bc03/csp_galaxev
source ./bc03.rm

# Assign environment
source csp_fits.tmp

# Build fits file
$bc03/buildfitsfile $name0 $ioptn
