"""process selection tkinter window."""

import tkinter


class ProcessSelector(object):
    """Process Selection class."""

    def __init__(self, master, list_like):
        """Initialize process selector class."""
        self.master = master
        self.v1 = tkinter.StringVar()
        self.lister(list_like)
        tkinter.Button(master, text="OK", command=self.ok_button).pack()

    def lister(self, items):
        """Create radio buttons from a list."""
        try:
            for text, mode in items:
                b = tkinter.Radiobutton(master, text=text,
                                        variable=self.v1, value=mode)
                b.select()
                b.pack()

        except Exception:
            for item in items:
                b = tkinter.Radiobutton(master, text=item, variable=self.v1,
                                        value=item)
                b.select()
                b.pack()

    def ok_button(self):
        """Destroy window when ok button is clicked."""
        self.master.destroy()


if __name__ == "__main__":

    # tuple_list = [("Thing one", "one"), ("Thing two", "two")]
    just_list = ["thing 1", "thing 2", "thing 3"]
    master = tkinter.Tk()
    # w = Window(master, tuple_list)
    w = ProcessSelector(master, just_list)
    master.mainloop()
    print(w.v1.get())
