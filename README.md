# PACESS

## Introduction

This is the repository of code for the PACESS paper. This contains a series of scripts which can be used for the pipeline described in the paper. The order of use of these scripts is shown in the diagram shown below.   

## Steps

```mermaid
flowchart TB

subgraph X[Sampling]
direction TB
    A[Tissue sample] --> B[Imaging]
end

subgraph Y[Cell-Extraction]
direction TB
    C[2D object detection neural network] --> H["Sort/arrange Data (02)"] --> D["3D prediction (03)"]
end

subgraph Z["Spatial-Analysis (05)\n"]
direction TB
    E[Spatial autocorrelation] --- F[Spatially-weighted regression model]
end

X --> Y -->I["Generate spatial dataframe (04)"] --> Z

X--->G["Extract boundaries (01)"] -.-> D["3D prediction (03)"]


classDef myscripts fill:#BDEDF2,stroke:#333,stroke-width:4px
classDef mysubscripts fill:#BDF2D9,stroke:#333
class G,H,D,I,Z myscripts
class E,F mysubscripts

```