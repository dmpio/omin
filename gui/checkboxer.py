"""process selection tkinter window."""

import tkinter


class CheckBoxer(object):
    """Process Selection class."""

    def __init__(self, master, list_like):
        """Initialize process selector class."""
        self.master = master
        self.item_dict = dict()
        self.lister(list_like)
        tkinter.Button(master, text="OK", command=self.ok_button).pack()

    def lister(self, items):
        """Create radio buttons from a list."""
        # self.item_dict = dict()
        for item in items:
            self.item_dict[item] = tkinter.Variable()
            b = tkinter.Checkbutton(master,
                                    text=item,
                                    variable=self.item_dict[item])
            # b.select()
            b.pack()

    def ok_button(self):
        """Destroy window when ok button is clicked."""
        self.master.destroy()


if __name__ == "__main__":

    # tuple_list = [("Thing one", "one"), ("Thing two", "two")]
    just_list = ["thing 1", "thing 2", "thing 3", "thing 4"]
    master = tkinter.Tk()
    # w = Window(master, tuple_list)
    w = CheckBoxer(master, just_list)
    master.mainloop()
