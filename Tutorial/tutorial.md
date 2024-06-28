# Tutorial

Here, I explain in details all the steps on this little program.

## Menu
When you execute `main.py`, a first window will be open to inform you to choose a working folder. Click on OK

![1st window](Images_tutorial/1.png)

Then, here you can choose a working directory or create one by clicking on the icon in top right (red arrow), enter a name and type enter on your keyboard and then select this new folder and click on **Open**.

![Menu](Images_tutorial/2.png)

After, the Menu window appears and here youi can choose 1 or all tasks and click on **Select**.

![Choose process directory](Images_tutorial/3.png)

## Smart-seq2 format

### Generate matrix
If you have check Smart-seq2 format, an information window will open to explain you to Choose folder(s) (_counts) to create a matrix if you have one file per cell, otherwise if you already have matrix, select nothing and click on Cancel on the next window and the processing will continue. This create a matrix with the number and the name of the counts folder in `Data`.
**Warning** your matrix have to have genes names in rows and barcodes/cell ID in columns. If you have genes IDs, you can use `modify_ids_names.py` and enter at the beginning of the script the name of your input file and the output name like this: **matrix_Number_Author.tab**.

![Window 2](Images_tutorial/4.png) ![Window 3](Images_tutorial/6.png)![Window 4](Images_tutorial/7.png)

### Create objects
The next step is creating object from the matrix. An information window will open and explain the nexte step.

![Window 4](Images_tutorial/8.png)

After, you select 1 or more matrix to create objects to pourchase processing and analysis.

![Window 5](Images_tutorial/10.png)

Object files will be in the working directory select or create at the begining in `Objects/` with the name `object_X_xxx_ori.h5ad`. It also generate a violin plot in `working_dir/Plots/X_xxx/` to see quality metrics of the dataset.

![Window 6](Images_tutorial/11.png)

  > [!NOTE] 
  > You can read .h5ad file with Seurat if you use also this packages in R.


