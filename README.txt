Usage:
./survival_interaction.R variables_input pivot_expression clinical_outcome output

variables_input:
Matrix of input variables, with columns as samples and rows as genes. This matrix can contain any type of information, including gene expression, mutation, and protein level.

Pivot_expression:
Matrix of pivot gene expression, with columns as samples and rows as genes. We set pivot genes as CD8A, CD8B, GZMA, GZMB, and PRF1. Please make sure that these genes should be profiled in the pivot expression matrix.

Clinical_outcome:
Response variables of survival and other clinical information. The first two columns should always be survival length and event status (death=1, alive=0). Other clinical information (e.g., age, gender, stage) could be included as later columns. The survival length could be either overall survival (OS) or progression-free survival (PFS).

Output:
There will be several types of output files with different postfix. Each file contains the z-score, p-value, and false discovery rate (FDR) for each gene tested.
      *.interaction: interaction scores between pivot (CTL level) and other genes with survival as the response.
      *.base:	     basic associations between the gene variable and survival outcome
      *.partial:     associations between the gene variable and survival outcome with the CTL level corrected in Cox-PH regression.
      *.main:	     the main effect associated with each gene in the interaction test.

Examples:
In the current folder, try the following commands.
   1: ./src/survival_interaction.R ./example/VanAllen.expression ./example/VanAllen.expression ./example/VanAllen.OS ./example/output_OS_Expression
   2: ./src/survival_interaction.R ./example/VanAllen.expression ./example/VanAllen.expression ./example/VanAllen.PFS ./example/output_PFS_Expression
   3: ./src/survival_interaction.R ./example/VanAllen.Mutation ./example/VanAllen.expression ./example/VanAllen.OS ./example/output_OS_Mutation
   4: ./src/survival_interaction.R ./example/VanAllen.Mutation ./example/VanAllen.expression ./example/VanAllen.PFS ./example/output_PFS_Mutation

Please contact: Peng Jiang (peng.jiang.software@gmail.com) if you have any questions or find some problems.
