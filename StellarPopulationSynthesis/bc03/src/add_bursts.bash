#!/bin/sh

# Script add_bursts.sh
source $bc03/.bc_bash
$bc03/add_bursts
source ./bc03.rm

# Assign environment
source csp_fits.tmp

# Build fits file
$bc03/buildfitsfile $name0 $ioptn
