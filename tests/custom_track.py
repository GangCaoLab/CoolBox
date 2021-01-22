from coolbox.core.track.base import Track


class CustomTrack(Track):
    def __init__(self):
        super().__init__(properties_dict={})

    def fetch_data(self, gr, **kwargs):
        return "Test Custom"

    def plot(self, ax, gr, **kwargs):
        x = gr.start + gr.length * 0.33
        ax.text(x, 0, self.fetch_data(gr), fontsize=50)
        ax.set_xlim(gr.start, gr.end)
