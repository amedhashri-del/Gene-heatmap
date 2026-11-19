import csv
import math
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from pathlib import Path


# This file opens a small desktop app for gene co-expression analysis.


def read_expression_file(file_name):
    genes = []
    values = []

    with open(file_name, "r", newline="") as file:
        reader = csv.reader(file)
        header = next(reader, None)

        if header is None or len(header) < 3:
            raise ValueError("The CSV needs a Gene column and at least two sample columns.")

        for row in reader:
            if not row:
                continue

            if len(row) != len(header):
                raise ValueError("One of the rows has a different number of columns.")

            gene_name = row[0].strip()

            if gene_name == "":
                raise ValueError("A gene name is missing.")

            gene_values = []

            for item in row[1:]:
                try:
                    gene_values.append(float(item))
                except ValueError:
                    raise ValueError("All expression values must be numbers.")

            genes.append(gene_name)
            values.append(gene_values)

    if len(genes) < 2:
        raise ValueError("Please include at least two genes.")

    if len(values[0]) < 2:
        raise ValueError("Please include at least two samples.")

    return genes, values


def pearson_correlation(first, second):
    first_mean = sum(first) / len(first)
    second_mean = sum(second) / len(second)

    top = 0
    first_bottom = 0
    second_bottom = 0

    for i in range(len(first)):
        first_difference = first[i] - first_mean
        second_difference = second[i] - second_mean

        top = top + first_difference * second_difference
        first_bottom = first_bottom + first_difference ** 2
        second_bottom = second_bottom + second_difference ** 2

    bottom = math.sqrt(first_bottom * second_bottom)

    if bottom == 0:
        return None

    return top / bottom


def make_correlation_matrix(values):
    matrix = []

    for first_gene in values:
        row = []

        for second_gene in values:
            result = pearson_correlation(first_gene, second_gene)
            row.append(result)

        matrix.append(row)

    return matrix


def correlation_colour(value):
    if value is None:
        return "#d9d9d9"

    value = max(-1, min(1, value))

    if value >= 0:
        red = 255
        green = int(255 * (1 - value))
        blue = int(255 * (1 - value))
    else:
        red = int(255 * (1 + value))
        green = int(255 * (1 + value))
        blue = 255

    return f"#{red:02x}{green:02x}{blue:02x}"


class GeneHeatmapApp:
    def __init__(self, root):
        self.root = root
        self.root.title("Gene Co-expression Heatmap")
        self.root.geometry("900x700")
        self.root.minsize(650, 500)

        self.genes = []
        self.values = []
        self.matrix = []
        self.file_name = None

        heading = tk.Label(
            root,
            text="Gene Co-expression Heatmap",
            font=("Arial", 20, "bold")
        )
        heading.pack(pady=15)

        description = tk.Label(
            root,
            text="Load a gene expression CSV file to calculate Pearson correlation.",
            font=("Arial", 11)
        )
        description.pack(pady=3)

        button_frame = tk.Frame(root)
        button_frame.pack(pady=12)

        tk.Button(
            button_frame,
            text="Load CSV File",
            command=self.load_file,
            width=16
        ).pack(side="left", padx=5)

        tk.Button(
            button_frame,
            text="Calculate Heatmap",
            command=self.calculate,
            width=18
        ).pack(side="left", padx=5)

        tk.Button(
            button_frame,
            text="Save Correlation CSV",
            command=self.save_matrix,
            width=20
        ).pack(side="left", padx=5)

        self.file_label = tk.Label(root, text="No file loaded", fg="gray")
        self.file_label.pack(pady=4)

        self.status_label = tk.Label(root, text="Load a CSV file to begin.", fg="blue")
        self.status_label.pack(pady=4)

        self.heatmap_frame = tk.Frame(root)
        self.heatmap_frame.pack(fill="both", expand=True, padx=15, pady=10)

        self.canvas = tk.Canvas(
            self.heatmap_frame,
            background="white",
            highlightthickness=1,
            highlightbackground="#cccccc"
        )
        self.canvas.pack(fill="both", expand=True)

        legend = tk.Label(
            root,
            text="Red = positive correlation     White = near zero     Blue = negative correlation",
            font=("Arial", 10)
        )
        legend.pack(pady=10)

        sample_file = Path(__file__).with_name("sample_data.csv")
        if sample_file.exists():
            self.load_data(str(sample_file))

    def load_data(self, file_name):
        try:
            self.genes, self.values = read_expression_file(file_name)
            self.file_name = file_name
            self.matrix = []
            self.file_label.config(text="Loaded: " + Path(file_name).name, fg="green")
            self.status_label.config(
                text=f"Loaded {len(self.genes)} genes and {len(self.values[0])} samples.",
                fg="blue"
            )
            self.canvas.delete("all")
        except Exception as error:
            messagebox.showerror("Could not load file", str(error))

    def load_file(self):
        file_name = filedialog.askopenfilename(
            title="Choose gene expression CSV",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )

        if file_name:
            self.load_data(file_name)

    def calculate(self):
        if not self.values:
            messagebox.showwarning("No data", "Please load a CSV file first.")
            return

        self.matrix = make_correlation_matrix(self.values)
        self.draw_heatmap()
        self.status_label.config(text="Heatmap calculated successfully.", fg="green")

    def draw_heatmap(self):
        self.canvas.delete("all")

        cell_size = 75
        left_space = 125
        top_space = 90

        width = left_space + len(self.genes) * cell_size + 20
        height = top_space + len(self.genes) * cell_size + 20

        self.canvas.config(scrollregion=(0, 0, width, height))

        for column, gene in enumerate(self.genes):
            x = left_space + column * cell_size + cell_size / 2
            self.canvas.create_text(
                x, top_space - 20,
                text=gene,
                font=("Arial", 9, "bold"),
                angle=45
            )

        for row, gene in enumerate(self.genes):
            y = top_space + row * cell_size + cell_size / 2
            self.canvas.create_text(
                left_space - 12, y,
                text=gene,
                font=("Arial", 9, "bold"),
                anchor="e"
            )

            for column, value in enumerate(self.matrix[row]):
                x1 = left_space + column * cell_size
                y1 = top_space + row * cell_size
                x2 = x1 + cell_size
                y2 = y1 + cell_size

                colour = correlation_colour(value)
                self.canvas.create_rectangle(
                    x1, y1, x2, y2,
                    fill=colour,
                    outline="white"
                )

                if value is None:
                    label = "N/A"
                else:
                    label = f"{value:.2f}"

                self.canvas.create_text(
                    (x1 + x2) / 2,
                    (y1 + y2) / 2,
                    text=label,
                    font=("Arial", 10)
                )

    def save_matrix(self):
        if not self.matrix:
            messagebox.showwarning("No result", "Calculate the heatmap first.")
            return

        save_name = filedialog.asksaveasfilename(
            title="Save correlation matrix",
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv")]
        )

        if not save_name:
            return

        with open(save_name, "w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Gene"] + self.genes)

            for i, gene in enumerate(self.genes):
                row = [gene]

                for value in self.matrix[i]:
                    if value is None:
                        row.append("N/A")
                    else:
                        row.append(round(value, 4))

                writer.writerow(row)

        messagebox.showinfo("Saved", "Correlation matrix saved successfully.")


root = tk.Tk()
app = GeneHeatmapApp(root)
root.mainloop()