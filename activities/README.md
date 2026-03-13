# Genome Segmentation Clustering Analysis
### CS4150 Project — Patrick McConnell

---

## 1. Project Overview

This project explores structural patterns in genomic data using unsupervised clustering techniques. The analysis focuses on **Genome Architecture Mapping (GAM)** segmentation data from the **Hist1 region on mouse chromosome 13**.

The goal is to cluster **nuclear profiles (NPs)** based on similarity in their genomic window presence and identify patterns of chromatin organization.

---

## 2. Motivation

Genome Architecture Mapping captures which genomic regions appear within individual nuclear slices. Comparing these slices across many nuclear profiles can reveal patterns of **3D genome organization**.

This project investigates whether nuclear profiles cluster into meaningful structural groups and whether those clusters correspond to specific chromatin features.

---

## 3. Data Processing

The dataset contains genomic windows and their presence across many nuclear profiles.

Each row represents a genomic window, and each column represents a nuclear profile.

Preprocessing steps include:

- Filtering the dataset to the **Hist1 region**
- Removing nuclear profiles with no detected windows
- Constructing a binary matrix representing genomic window presence

---

## 4. Similarity and Clustering

Similarity between nuclear profiles is computed using the **Jaccard index**, which measures the overlap between binary vectors.

A distance matrix derived from this similarity is used to perform **k-medoids clustering**. Cluster quality is evaluated using **within-cluster variation**.

---

## 5. Visualization

Clusters are visualized using **heatmaps** of the segmentation matrix.

- Rows represent nuclear profiles  
- Columns represent genomic windows  
- Cells represent presence or absence of a genomic region  

These visualizations help reveal structural patterns in the genome.

---

## 6. Tools

The project is implemented in **Python** using:

- pandas  
- numpy  
- matplotlib  
- seaborn  

---

## 7. Running the Project

To run the analysis pipeline:
- python main.py

The script loads the dataset, computes similarity matrices, performs clustering, and generates visualizations.
