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
If you have check Smart-seq2 format, an information window will open to explain you to Choose folder(s) (_counts) to create a matrix if you have one file per cell, otherwise if you already have matrix, select nothing and click on Cancel on the next window and the processing will continue. **Warning** your matrix have to have genes names in rows and barcodes/cell ID in columns. If you have genes IDs, you can use `modify_ids_names.py` and enter at the beginning of the script the name of your input file and the output name like this: **matrix_Number_Author.tab**.

![Window 2](Images_tutorial/4.png) ![Window 3](Images_tutorial/6.png)


