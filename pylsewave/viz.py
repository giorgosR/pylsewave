"""
pylsewave module with vizualisation utils.
"""
import matplotlib.pyplot as plt
from matplotlib import rcParams
import numpy as np
import glob
import os
import shutil
import time

__author__ = "Georgios E. Ragkousis"
__notes__ = "Class body recycled from Hans Peter Langtangen  (code not used in any example)"
rcParams['font.family'] = 'serif'
rcParams['font.size'] = 16

WRITE_ARRAY = 0


def viz(I, V, f, U_0, U_L, c, L, dt, C, T, umin, umax, solver, animate=True,
        tool="matplotlib"):
    def plot_u_st(u, x, t, n):
        """user_action function for solver."""
        plt.plot(x, u, 'r-',
                 xlabel='x', ylabel='u',
                 axis=[0, L, umin, umax],
                 title='t=%f' % t[n], show=True)
        # Let the initial condition stay on the screen for 2
        # seconds, else insert a pause of 0.2 s between each plot
        time.sleep(2) if t[n] == 0 else time.sleep(0.2)
        plt.savefig('frame_%04d.png' % n)  # for movie making

    class PlotMatplotlib:
        def __call__(self, u, x, t, n):
            """user_action function for solver."""
            if n == 0:
                plt.ion()
                self.lines = plt.plot(x, u, 'r-')
                plt.xlabel('x')
                plt.ylabel('u')
                plt.axis([0, L, umin, umax])
                plt.legend(['t=%f' % t[n]], loc='lower left')
            else:
                self.lines[0].set_ydata(u)
                plt.legend(['t=%f' % t[n]], loc='lower left')
                plt.draw()
            time.sleep(2) if t[n] == 0 else time.sleep(0.2)
            plt.savefig('tmp_%04d.png' % n)  # for movie making

    if tool == 'matplotlib':
        import matplotlib.pyplot as plt
        plot_u = PlotMatplotlib()
    elif tool == 'scitools':
        import scitools.std as plt  # scitools.easyviz interface
        plot_u = plot_u_st
    import time, glob, os

    # Clean up old movie frames
    for filename in glob.glob('tmp_*.png'):
        os.remove(filename)

    # Call solver and do the simulation
    user_action = plot_u if animate else None
    u, x, t, cpu_ = solver(I, V, f, U_0, U_L, c, L, dt, C, T, user_action)

    # Make video files
    fps = 4  # frames per second
    codec2ext = dict(flv='flv', libx264='mp4', libvpx='webm',
                     libtheora='ogg')  # video formats
    filespec = 'tmp_%04d.png'
    movie_program = 'ffmpeg'  # or 'avconv'
    for codec in codec2ext:
        ext = codec2ext[codec]
        cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s ' \
              '-vcodec %(codec)s movie.%(ext)s' % vars()
        os.system(cmd)

    if tool == 'scitools':
        # Make an HTML play for showing the animation in a browser
        plt.movie('tmp_*.png', encoder='html', fps=fps,
                  output_file='movie.html')
    return cpu_


def merge_zip_archives(individual_archives, archive_name):
    """
    Merge individual zip archives made with numpy.savez into
    one archive with name archive_name.
    The individual archives can be given as a list of names
    or as a Unix wild chard filename expression for glob.glob.
    The result of this function is that all the individual
    archives are deleted and the new single archive made.
    """
    import zipfile
    with zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED,
                         allowZip64=True) as archive:
        #    archive = zipfile.ZipFile(archive_name, 'w', zipfile.ZIP_DEFLATED,
        #                              allowZip64=True)
        if isinstance(individual_archives, (list, tuple)):
            filenames = individual_archives
        elif isinstance(individual_archives, str):
            filenames = glob.glob(individual_archives)

        # Open each archive and write to the common archive
        for filename in filenames:
            with zipfile.ZipFile(filename, 'r', zipfile.ZIP_DEFLATED) as f:
                #        f = zipfile.ZipFile(filename,  'r',
                #                            zipfile.ZIP_DEFLATED)
                for name in f.namelist():
                    data = f.open(name, 'r')
                    # Save under name without .npy
                    archive.writestr(name[:-4], data.read())
                f.close()
                #            try:
                #                os.remove(filename)
                #            except WindowsError:
                #                pass
        archive.close()
    print("[-] Cleaning temporary files...")
    for filename in filenames:
        try:
            os.remove(filename)
        except WindowsError:
            pass


class PlotAndStoreSolution:
    """
    Class for the user_action function in solver.
    Visualizes the solution only.
    """

    def __init__(
            self,
            casename='tmp',  # Prefix in filenames
            umin=-1, umax=1,  # Fixed range of y axis
            pause_between_frames=None,  # Movie speed
            backend='matplotlib',  # or 'gnuplot' or None
            screen_movie=True,  # Show movie on screen?
            title='',  # Extra message in title
            skip_frame=1,  # Skip every skip_frame frame
            filename=None):  # Name of file with solutions
        self.casename = casename
        self.yaxis = [umin, umax]
        self.pause = pause_between_frames
        self.backend = backend
        if backend is None:
            # Use native matplotlib
            import matplotlib.pyplot as plt
        elif backend in ('matplotlib', 'gnuplot'):
            module = 'scitools.easyviz.' + backend + '_'
            exec ('import %s as plt' % module)
        self.plt = plt
        self.screen_movie = screen_movie
        self.title = title
        self.skip_frame = skip_frame
        self.filename = filename
        self.index = 0
        if filename is not None:
            # Store time points when u is written to file
            self.t = []
            filenames = glob.glob('.' + self.filename + '*.dat.npz')
            for filename in filenames:
                os.remove(filename)

        # Clean up old movie frames
        for filename in glob.glob('frame_*.png'):
            os.remove(filename)

    def __call__(self, u, x, t, n, verbose=True, store_data=False, vessel_index=None):
        """
        Callback function user_action, call by solver:
        Store solution, plot on screen and save to file.
        """
        if n % self.skip_frame != 0:
            return
        if verbose is True:
            np.set_printoptions(precision=3, suppress=True)
            print("Vessel %d"  %vessel_index + "u: ", u)
        if store_data is True:
            # Save solution u to a file using numpy.savez
            if self.filename is not None and vessel_index is None:
                name = 'u%010d' % n  # array name
                kwargs = {name: u}
                fname = '.' + self.filename + '_' + name + '.dat'
                # store solution to file
                np.savez(fname, **kwargs)
                self.t.append(t[n])  # store corresponding time value
                if n == 0:  # save x once
                    np.savez('.' + self.filename + '_x.dat', x=x)

            elif self.filename is not None and vessel_index is not None:
                name = 'u%010d_vessel_%03d' % (n, vessel_index)  # array name
                kwargs = {name: u}
                fname = '.' + self.filename + '_' + name + '.dat'
                # store solution to file
                np.savez(fname, **kwargs)
                if vessel_index == 0:
                    self.t.append(t[n])  # store corresponding time value
                if n == 0:  # save x once
                    namex = 'x_vessel_%03d' % (vessel_index)
                    kwargs = {namex: x}
                    fnamex = '.' + self.filename + '_' + namex + '.dat'
                    np.savez(fnamex, **kwargs)

            # # Animate
            # # if n % self.skip_frame != 0:
            # #    return
            # title = 't=%.3f' % t[n]
            # if self.title:
            #     title = self.title + ' ' + title
            # if self.backend is None:
            #     # native matplotlib animation
            #     if n == 0:
            #         self.plt.ion()
            #         self.lines = self.plt.plot(x, u[WRITE_ARRAY], 'r-')
            #         self.plt.axis([x[0], x[-1],
            #                        self.yaxis[0], self.yaxis[1]])
            #         self.plt.xlabel('x')
            #         self.plt.ylabel('u')
            #         self.plt.title(title)
            #         self.plt.legend(['t=%.3f' % t[n]])
            #     else:
            #         # Update new solution
            #         self.lines[0].set_ydata(u[WRITE_ARRAY])
            #         self.plt.legend(['t=%.3f' % t[n]])
            #         self.plt.draw()
            # else:
            #     # scitools.easyviz animation
            #     self.plt.plot(x, u, 'r-',
            #                   xlabel='x', ylabel='u',
            #                   axis=[x[0], x[-1],
            #                         self.yaxis[0], self.yaxis[1]],
            #                   title=title,
            #                   show=self.screen_movie)
            # # pause
            # if t[n] == 0:
            #     time.sleep(2)  # let initial condition stay 2 s
            # else:
            #     if self.pause is None:
            #         pause = 0.2 if u[WRITE_ARRAY].size < 100 else 0
            #     time.sleep(pause)
            # if vessel_index is None:
            #     self.plt.savefig('frame_%03d.png' % (n/self.skip_frame))
            #     # this is to enumerate the frames
            #     self.index += 1
            # else:
            #     self.plt.savefig('frame_%03d_vessel_%03d.png' % (n/self.skip_frame, vessel_index ))

    def make_movie_file(self):
        """
        Create subdirectory based on casename, move all plot
        frame files to this directory, and generate
        an index.html for viewing the movie in a browser
        (as a sequence of PNG files).
        """
        # Make HTML movie in a subdirectory
        directory = self.casename

        if os.path.isdir(directory):
            shutil.rmtree(directory)  # rm -rf directory
        os.mkdir(directory)  # mkdir directory
        # mv frame_*.png directory
        for filename in glob.glob('frame_*.png'):
            os.rename(filename, os.path.join(directory, filename))
        os.chdir(directory)  # cd directory

        fps = 24  # frames per second
        if self.backend is not None:
            from scitools.std import movie
            movie('frame_*.png', encoder='html',
                  output_file='index.html', fps=fps)

        # Make other movie formats: Flash, Webm, Ogg, MP4
        codec2ext = dict(flv='flv', libx264='mp4', libvpx='webm',
                         libtheora='ogg')
        filespec = 'frame_%03d.png'
        movie_program = 'ffmpeg'  # or 'avconv'
        for codec in codec2ext:
            ext = codec2ext[codec]
            cmd = '%(movie_program)s -r %(fps)d -i %(filespec)s ' \
                  '-vcodec %(codec)s movie.%(ext)s' % vars()
            os.system(cmd)

        os.chdir(os.pardir)  # move back to parent directory

    def close_file(self, hashed_input):
        """
        Merge all files from savez calls into one archive.
        hashed_input is a string reflecting input data
        for this simulation (made by solver).
        """
        if self.filename is not None:
            # Save all the time points where solutions are saved
            np.savez('.' + self.filename + '_t.dat',
                     t=np.array(self.t, dtype=float))

            # Merge all savez files to one zip archive
            archive_name = '.' + hashed_input + '_archive.npz'
            filenames = glob.glob('.' + self.filename + '*.dat.npz')
            merge_zip_archives(filenames, archive_name)
            print('Archive name:', archive_name)
            # data = numpy.load(archive); data.files holds names
            # data[name] extract the array


class PlotMediumAndSolution(PlotAndStoreSolution):
    def __init__(self, medium, **kwargs):
        """Mark medium in plot: medium=[x_L, x_R]."""
        self.medium = medium
        PlotAndStoreSolution.__init__(self, **kwargs)

    def __call__(self, u, x, t, n):
        # Save solution u to a file using numpy.savez
        if self.filename is not None:
            name = 'u%04d' % n  # array name
            kwargs = {name: u}
            fname = '.' + self.filename + '_' + name + '.dat'
            np.savez(fname, **kwargs)
            self.t.append(t[n])  # store corresponding time value
            if n == 0:  # save x once
                np.savez('.' + self.filename + '_x.dat', x=x)

        # Animate
        if n % self.skip_frame != 0:
            return
        # Plot u and mark medium x=x_L and x=x_R
        x_L, x_R = self.medium
        umin, umax = self.yaxis
        title = 'Nx=%d' % (x.size - 1)
        if self.title:
            title = self.title + ' ' + title
        if self.backend is None:
            # native matplotlib animation
            if n == 0:
                self.plt.ion()
                self.lines = self.plt.plot(
                    x, u, 'r-',
                    [x_L, x_L], [umin, umax], 'k--',
                    [x_R, x_R], [umin, umax], 'k--')
                self.plt.axis([x[0], x[-1],
                               self.yaxis[0], self.yaxis[1]])
                self.plt.xlabel('x')
                self.plt.ylabel('u')
                self.plt.title(title)
                self.plt.text(0.75, 1.0, 'c lower')
                self.plt.text(0.32, 1.0, 'c=1')
                self.plt.legend(['t=%.3f' % t[n]])
            else:
                # Update new solution
                self.lines[0].set_ydata(u)
                self.plt.legend(['t=%.3f' % t[n]])
                self.plt.draw()
        else:
            # scitools.easyviz animation
            self.plt.plot(x, u, 'r-',
                          [x_L, x_L], [umin, umax], 'k--',
                          [x_R, x_R], [umin, umax], 'k--',
                          xlabel='x', ylabel='u',
                          axis=[x[0], x[-1],
                                self.yaxis[0], self.yaxis[1]],
                          title=title,
                          show=self.screen_movie)
        # pause
        if t[n] == 0:
            time.sleep(2)  # let initial condition stay 2 s
        else:
            if self.pause is None:
                pause = 0.2 if u.size < 100 else 0
            time.sleep(pause)

        self.plt.savefig('frame_%04d.png' % (n))

        if n == (len(t) - 1):  # finished with this run, close plot
            self.plt.close()

if __name__ == '__main__':
    raise NotImplementedError('Module is not idented for direct execution')