from scipy.sparse import csr_matrix

if __name__ == "__main__":
    import tkinter as tk
    import subprocess
    import anndata
    import cell_extract
    import compare_expr
    import farell_new
    import cross_reference
    import pandas as pd
    import scanpy as sc

    def search_gene():
        adata = sc.read_mtx("C:\\Users\\rohan\OneDrive\Documents\gene_analysis\\matrix.mtx").T

        adata.obs = pd.read_csv("C:\\Users\\rohan\OneDrive\Documents\gene_analysis\\barcodes.tsv",
                                index_col=0, sep="\t")
        adata.var = pd.read_csv("C:\\Users\\rohan\OneDrive\Documents\gene_analysis\\features.tsv",
                                index_col=1,
                                sep="\t")
        print("loaded")
        # print("start")
        # adata = adata[adata.obs['timepoint'] == selected_option.get()]
        # print("Done")
        # Clear the previous result

        result_text.delete('1.0', tk.END)

        # Get the gene name from the search bar

        gene_names = search_entry.get().replace(" ", "").split(",")
        gene_levels = []
        count = 0
        for i in gene_names:
            gene_names[count] = i[1:len(i)]
            count += 1
            gene_levels.append(i[0])

        # utilize cell extract
        print("start")
        internal, external = cell_extract.extract(adata, gene_names, gene_levels)


        # use farell_new
        loc1 = farell_new.create_genes(adata, internal)
        # use cross_reference
        loc2 = cross_reference.create_gene(adata, external)

        # use compare_expr
        s1 = max(loc1.values())
        s2 = max(loc2.values())
        top_100 = compare_expr.results(loc1, loc2, s1, s2)
        bottom_100 = top_100
        bottom_100.reverse()


        # # Process the gene name and store the result in a variable (replace this with your own logic)
        # with open('results.txt', 'r') as file:
        #     data = file.readlines()
        result = f"Results for {gene_names}:\n"
        t = len(top_100)
        for i in range(t):
            result += str(
                " " * 15 + top_100[i].replace("\n", "") + (30 - len(top_100[i])) * " " + bottom_100[i].replace("\n",
                                                                                                               "") + "\n")
        # Display the result in the text area
        result_text.insert(tk.END, result)


    # Create the tkinter app
    app = tk.Tk()

    # # Load the background image
    # background_image = tk.PhotoImage(file="DNA_background.png")
    # #
    # # # Set the background image on a label
    # background_label = tk.Label(app, image=background_image)
    # background_label.place(x=0, y=0, relwidth=1, relheight=1)
    selected_option = tk.StringVar(app)
    selected_option.set('10hpf')  # Set the default selected option

    # Create the dropdown menu
    dropdown_menu = tk.OptionMenu(app, selected_option, '10hpf', '12hpf', '14hpf', '16hpf', '19hpf', '1dpf', '2dpf',
                                  '3dpf', '5dpf', '10dpf', "full")
    dropdown_menu.pack()
    # Create the search bar
    search_entry = tk.Entry(app)
    search_entry.pack()
    title_label = tk.Label(app, text="Gene Search Tool - RD", font=("Helvetica", 18, "bold"), fg="black")
    title_label.pack(pady=10)
    # Create the search button
    search_button = tk.Button(app, text="Search", command=lambda: search_gene())
    search_button.pack()

    # Create the text area to display the results
    result_text = tk.Text(app)
    result_text.pack()


    def update_time():
        global adata
        timepoint = selected_option.get()
        adata = adata[adata.obs['timepoint'] == timepoint]


    # Create a new button here
    update_button = tk.Button(app, text="Update Time", command=update_time)
    update_button.pack(side=tk.LEFT, padx=5)
    # Create the additional text
    additional_text = tk.Label(app, text="Click Download for additional data")
    additional_text.pack()

    # Run the gene_similarity.py script when the app starts

    # Start the tkinter event loop
    app.mainloop()
