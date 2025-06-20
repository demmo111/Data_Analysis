# 🧬 Unraveling the Genetic and Pathway Overlap Between SARS-CoV-2 and Alzheimer’s Disease

This repository contains the complete analysis and results from a data-driven bioinformatics research project aimed at uncovering shared genetic markers and biological pathways between **COVID-19 (SARS-CoV-2)** and **Alzheimer’s Disease (AD)**.

---

## 📌 Project Title

**"Unraveling the Genetic and Pathway Overlap Between SARS-CoV-2 and Alzheimer’s Disease: Insights for Future Therapeutic Strategies"**

---

## 👥 Contributors

- Rana Mohamed Kamel — 202001118  
- Reem Sharaf Eldeen Hassan — 211001887  
- Mohamed Abdelkader Ragab — 221001955  
- Alaa Ashraf Lotfy — 221001804  
- Supervised by: Dr. Mohamed El Sayeh

---

## 🧪 Objective

To explore the genetic intersection and overlapping molecular mechanisms between COVID-19 and Alzheimer's Disease using bioinformatics tools, in order to uncover shared therapeutic targets and better understand the impact of viral infections on neurodegeneration.

---

## 🔬 Methods

1. **Data Collection:**  
   - Publicly available RNA-seq datasets from **GEO**:  
     - COVID-19: GSE171110  
     - Alzheimer's Disease: GSE53480  
   - Selection of 10 disease and 10 control samples from each dataset.

2. **Differential Expression Analysis (DEA):**  
   - Performed using **GEO2R** to obtain DEGs.

3. **Common Gene Identification:**  
   - DEGs from both datasets intersected using R.

4. **Functional Annotation and Enrichment Analysis:**  
   - Tools: **ShinyGO**, **KEGG** database.
   - Analyzed roles in inflammation, immune response, mitochondrial dysfunction.

5. **Validation:**  
   - Applied adjusted p-values and FDR to ensure statistical reliability.

---

## 🧬 Key Findings

- **Common genes identified:** C3, C2, ND5, F5, F12, F8
- These genes are involved in:  
  - Complement & coagulation cascades  
  - Mitochondrial dysfunction  
  - Immune dysregulation

- **Enriched KEGG pathways:**  
  - hsa04610: Complement and Coagulation Cascades  
  - hsa05171: Coronavirus Disease (COVID-19)  
  - hsa05012: Parkinson's Disease  
  - hsa05016: Huntington's Disease  
  - hsa05322: Systemic Lupus Erythematosus

---

## 📈 Results Highlights

- **Volcano plots**, **UMAP**, **Boxplots**, and **Venn Diagrams** were generated for visualization.  
- Clear **genetic overlap** between SARS-CoV-2 and AD identified.  
- Pathway analysis showed that **inflammation, coagulation, and mitochondrial dysfunction** are shared hallmarks.

---

## 💡 Implications

- Potential shared therapeutic targets were found (e.g., targeting complement activation or mitochondrial function).
- Suggests COVID-19 may accelerate Alzheimer’s progression.
- Proposes **dual-intervention strategies** for treatment.

---

## 📚 References

Key sources include datasets from [GEO](https://www.ncbi.nlm.nih.gov/geo/), tools like [ShinyGO](http://bioinformatics.sdstate.edu/go/), and literature from PubMed and Nature Scientific Reports.

---


---

## 📄 License

This project is developed as part of a university coursework and is intended for educational and research purposes.


