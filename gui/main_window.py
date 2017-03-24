import tkinter
from tkinter import filedialog

class omin_app(tkinter.Tk):
    def __init__(self):
        self.selected_files = None
        self.omin_obj = None

    def run(self):
        root = tkinter.Tk()
        root.title(('Omin GUI Demo'))
        root.resizable(width=True, height=True)
        root.minsize(width=600, height=200)
        # root.maxsize(width=666, height=666)
        root.call('wm', 'attributes', '.', '-topmost', True)
        tkinter.Button(root, text=('Select Input'), command=self.select_files, width=10).pack(side=tkinter.TOP)
        tkinter.Listbox(root).pack(side=tkinter.BOTTOM)
        root.mainloop()


    def select_files(self):
        selected_files = filedialog.askopenfilename(multiple=True)
        self.selected_files = selected_files
        return


if __name__ == "__main__":
    oa = omin_app()
    oa.run()
