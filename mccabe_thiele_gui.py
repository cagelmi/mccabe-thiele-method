import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
import json # For Save/Load parameters

class McCabeThieleApp:
    # Class constants for tolerances and limits
    CLS_TOLERANCE_EQ = 1e-6  # For equality checks with mole fractions
    CLS_TOLERANCE_DIV = 1e-9 # For preventing division by zero
    CLS_MAX_PLATES = 150     # Max plates safety break
    CLS_HOVER_DIST_FACTOR = 0.02 # Factor for hover sensitivity (fraction of axis width)


    def __init__(self, root):
        self.root = root
        self.root.title("MCT")
        # self.root.geometry("850x700") # Adjusted size for new elements

        # --- Store default parameters ---
        self.default_params = {
            "zf": "0.50", "q": "0.80", "xd": "0.95",
            "xb": "0.10", "R": "1.8", "alpha": "2.7"
        }

        # --- Initialize instance variables for results ---
        self.xi = np.nan
        self.yi = np.nan
        self.feed_plate_number = "N/A"
        self.current_num_plates = "N/A"

        # --- Menu ---
        menubar = tk.Menu(root)
        helpmenu = tk.Menu(menubar, tearoff=0)
        helpmenu.add_command(label="About", command=self.show_about)
        menubar.add_cascade(label="Help", menu=helpmenu)
        root.config(menu=menubar)

        # --- Main Title ---
        title_label = ttk.Label(root, text="McCabe and Thiele Graphical Method", font=("Arial", 16, "bold"))
        title_label.pack(pady=10)

        main_paned_window = ttk.PanedWindow(root, orient=tk.HORIZONTAL)
        main_paned_window.pack(fill=tk.BOTH, expand=True, padx=10, pady=5)

        plot_frame = ttk.Frame(main_paned_window, width=550, height=550)
        plot_frame.pack_propagate(False)
        main_paned_window.add(plot_frame, weight=2)

        self.fig, self.ax = plt.subplots()
        self.canvas = FigureCanvasTkAgg(self.fig, master=plot_frame)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill=tk.BOTH, expand=True)
        
        toolbar_frame = ttk.Frame(plot_frame)
        toolbar_frame.pack(fill=tk.X, side=tk.BOTTOM)
        toolbar = NavigationToolbar2Tk(self.canvas, toolbar_frame)
        toolbar.update()

        self.ax.set_xlabel("X", fontweight="bold")
        self.ax.set_ylabel("Y", fontweight="bold")
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        self.ax.plot([0, 1], [0, 1], 'k-', lw=1)
        self.ax.grid(True, linestyle=':', alpha=0.7)
        
        # --- Hover Information Setup ---
        self.hover_text = self.ax.text(0, 0, "", va="bottom", ha="left",
                                          bbox=dict(boxstyle="round,pad=0.3", fc="lemonchiffon", alpha=0.85),
                                          visible=False, zorder=10, fontsize=8)
        self.plot_hover_points = [] 
        self.canvas.mpl_connect('motion_notify_event', self._on_motion)
        self.canvas.draw()


        controls_frame_outer = ttk.Frame(main_paned_window, width=300)
        controls_frame_outer.pack_propagate(False)
        main_paned_window.add(controls_frame_outer, weight=1)

        controls_frame = ttk.Frame(controls_frame_outer)
        controls_frame.pack(padx=10, pady=10, fill=tk.BOTH, expand=True)

        op_conditions_frame = ttk.LabelFrame(controls_frame, text="Operation conditions")
        op_conditions_frame.pack(fill=tk.X, pady=5)

        self.inputs = {}
        conditions_spec = [
            ("Feed mole fraction (zf):", "zf"), ("Feed quality (q):", "q"),
            ("Distillate mole fraction (xd):", "xd"), ("Bottom mole fraction (xb):", "xb"),
            ("Reflux ratio (R):", "R")
        ]

        for i, (label_text, key) in enumerate(conditions_spec):
            ttk.Label(op_conditions_frame, text=label_text).grid(row=i, column=0, padx=5, pady=2, sticky="w")
            entry = ttk.Entry(op_conditions_frame, width=10)
            entry.grid(row=i, column=1, padx=5, pady=2, sticky="ew")
            entry.insert(0, self.default_params[key])
            self.inputs[key] = entry
        op_conditions_frame.columnconfigure(1, weight=1)

        rel_vol_frame = ttk.LabelFrame(controls_frame, text="Relative volatility")
        rel_vol_frame.pack(fill=tk.X, pady=5)
        ttk.Label(rel_vol_frame, text="alpha:").grid(row=0, column=0, padx=5, pady=2, sticky="w")
        entry_alpha = ttk.Entry(rel_vol_frame, width=10)
        entry_alpha.grid(row=0, column=1, padx=5, pady=2, sticky="ew")
        entry_alpha.insert(0, self.default_params["alpha"])
        self.inputs["alpha"] = entry_alpha
        rel_vol_frame.columnconfigure(1, weight=1)

        # --- Calculation Results Display ---
        results_frame = ttk.LabelFrame(controls_frame, text="Calculation Results")
        results_frame.pack(fill=tk.X, pady=10)

        self.num_plates_label = ttk.Label(results_frame, text="Number of plates: --", font=("Arial", 10))
        self.num_plates_label.pack(pady=2, anchor="w", padx=5)
        
        self.feed_plate_label = ttk.Label(results_frame, text="Feed plate number: --", font=("Arial", 10))
        self.feed_plate_label.pack(pady=2, anchor="w", padx=5)

        self.xi_label = ttk.Label(results_frame, text="Intersection xi: --", font=("Arial", 10))
        self.xi_label.pack(pady=2, anchor="w", padx=5)
        
        self.yi_label = ttk.Label(results_frame, text="Intersection yi: --", font=("Arial", 10))
        self.yi_label.pack(pady=2, anchor="w", padx=5)


        # --- Buttons Frame ---
        button_frame = ttk.Frame(controls_frame)
        button_frame.pack(pady=10, fill=tk.X)

        self.calc_button = ttk.Button(button_frame, text="Calculate", command=self.calculate_and_plot)
        self.calc_button.pack(side=tk.LEFT, padx=2, expand=True, fill=tk.X)

        self.export_button = ttk.Button(button_frame, text="Export Plot", command=self.export_plot)
        self.export_button.pack(side=tk.LEFT, padx=2, expand=True, fill=tk.X)
        
        # --- Parameter Management Buttons Frame ---
        param_button_frame = ttk.Frame(controls_frame)
        param_button_frame.pack(pady=5, fill=tk.X)

        self.load_params_button = ttk.Button(param_button_frame, text="Load Params", command=self._load_parameters)
        self.load_params_button.pack(side=tk.LEFT, padx=2, expand=True, fill=tk.X)
        
        self.save_params_button = ttk.Button(param_button_frame, text="Save Params", command=self._save_parameters)
        self.save_params_button.pack(side=tk.LEFT, padx=2, expand=True, fill=tk.X)

        self.reset_button = ttk.Button(param_button_frame, text="Reset Defaults", command=self._reset_to_defaults)
        self.reset_button.pack(side=tk.LEFT, padx=2, expand=True, fill=tk.X)
        
        author_label = ttk.Label(root, text="by Claudio A. Gelmi. (Python adaptation)")
        author_label.pack(side=tk.LEFT, padx=10, pady=5)

        self.calculate_and_plot(first_run=True)

    def _get_inputs(self):
        values = {}
        try:
            for key, entry_widget in self.inputs.items():
                values[key] = float(entry_widget.get())
            return values
        except ValueError:
            messagebox.showerror("Input Error", "All input fields must be valid numbers.")
            return None

    def _validate_inputs(self, params):
        zf, q, xd, xb, R, alpha = params['zf'], params['q'], params['xd'], params['xb'], params['R'], params['alpha']
        if not (0 < zf < 1 and 0 < xb < 1 and 0 < xd < 1):
            messagebox.showerror("Input Error", "Molar fractions (zf, xb, xd) must be between 0 and 1 (exclusive).")
            return False
        if not (xb < zf < xd):
             messagebox.showerror("Input Error", "Ensure xb < zf < xd for standard distillation.")
             return False
        if alpha <= 1:
            messagebox.showerror("Input Error", "Relative volatility (alpha) must be greater than 1.")
            return False
        if R <= 0:
            messagebox.showerror("Input Error", "Reflux ratio (R) must be positive.")
            return False
        return True

    def equilib_y_from_x(self, x, alpha):
        return alpha * x / (1 + (alpha - 1) * x)

    def equilib_x_from_y(self, y, alpha):
        denominator = alpha - y * (alpha - 1)
        if abs(denominator) < self.CLS_TOLERANCE_DIV:
            if y > 0.999: return 1.0
            if y < 0.001: return 0.0
            return np.nan
        return y / denominator
        
    def _reset_results_display(self):
        self.num_plates_label.config(text="Number of plates: --")
        self.feed_plate_label.config(text="Feed plate number: --")
        self.xi_label.config(text="Intersection xi: --")
        self.yi_label.config(text="Intersection yi: --")
        self.xi, self.yi = np.nan, np.nan
        self.feed_plate_number = "N/A"
        self.current_num_plates = "N/A"
        self.plot_hover_points.clear()


    def calculate_and_plot(self, first_run=False):
        self._reset_results_display() 

        params = self._get_inputs()
        if not params:
            self.ax.clear()
            self._setup_plot_labels_and_legend()
            self.canvas.draw()
            return

        if not self._validate_inputs(params) and not first_run:
            self.ax.clear()
            self._setup_plot_labels_and_legend()
            self.canvas.draw()
            return
        
        zf, q_val, xd, xb, R_val, alpha = params['zf'], params['q'], params['xd'], params['xb'], params['R'], params['alpha']

        self.ax.clear()
        self.plot_hover_points.clear()

        x_eq_pts = np.linspace(0, 1, 200)
        y_eq_pts = self.equilib_y_from_x(x_eq_pts, alpha)
        self.ax.plot(x_eq_pts, y_eq_pts, 'r-', lw=1.5, label="Equilibrium Curve")
        for x_pt, y_pt in zip(x_eq_pts[::10], y_eq_pts[::10]): 
             self.plot_hover_points.append({'x': x_pt, 'y': y_pt, 'label': f"Eq: ({x_pt:.3f}, {y_pt:.3f})", 'type': 'eq'})

        self.ax.plot([0, 1], [0, 1], 'k-', lw=1, label="y=x")

        if abs(q_val - 1.0) < self.CLS_TOLERANCE_EQ: 
            self.xi = zf
            self.yi = (R_val / (R_val + 1)) * self.xi + xd / (R_val + 1)
        else: 
            m_q = q_val / (q_val - 1)
            c_q = zf * (1 - m_q) 
            m_rol = R_val / (R_val + 1)
            c_rol = xd / (R_val + 1)
            
            if abs(m_q - m_rol) < self.CLS_TOLERANCE_DIV:
                 if not first_run: messagebox.showerror("Calculation Error", "Operating lines are parallel. Check inputs.")
                 self._reset_results_display()
                 self._setup_plot_labels_and_legend()
                 self.canvas.draw()
                 return
            
            self.xi = (c_rol - c_q) / (m_q - m_rol)
            self.yi = m_rol * self.xi + c_rol
        
        self.plot_hover_points.append({'x': self.xi, 'y': self.yi, 
                                       'label': f"Intersection (xi, yi)\nx={self.xi:.3f}, y={self.yi:.3f}", 
                                       'type': 'intersect'})

        y_eq_at_xi = self.equilib_y_from_x(self.xi, alpha)
        if self.yi > y_eq_at_xi + self.CLS_TOLERANCE_EQ :
            if not first_run: messagebox.showerror("Calculation Error", "Distillation not possible (Operating lines cross above equilibrium at feed). Try different R or other conditions.")
            self._reset_results_display()
            self._setup_plot_labels_and_legend()
            self.canvas.draw()
            return

        self.ax.plot([xd, self.xi], [xd, self.yi], 'b-', lw=1, label="ROL")
        if abs(self.xi - xb) < self.CLS_TOLERANCE_DIV:
            if not first_run: messagebox.showerror("Calculation Error", "xi is too close to xb, cannot determine SOL slope.")
            self._reset_results_display()
            self._setup_plot_labels_and_legend()
            self.canvas.draw()
            return
        ss = (self.yi - xb) / (self.xi - xb) 
        self.ax.plot([xb, self.xi], [xb, self.yi], 'b-', lw=1, label="SOL")
        self.ax.plot([zf, self.xi], [zf, self.yi], 'b--', lw=1, label="q-line")

        num_plates = 0
        self.feed_plate_number = "N/A"
        x_curr, y_curr = xd, xd
        
        self.plot_hover_points.append({'x': x_curr, 'y': y_curr, 'label': "Start (xd,xd)", 'type': 'stage'})

        while x_curr > xb + self.CLS_TOLERANCE_EQ and num_plates < self.CLS_MAX_PLATES :
            num_plates += 1
            
            x_on_equilibrium = self.equilib_x_from_y(y_curr, alpha)
            
            self.ax.plot([x_curr, x_on_equilibrium], [y_curr, y_curr], 'm-', lw=0.8)
            self.ax.text(x_on_equilibrium - 0.01, y_curr + 0.02, str(num_plates), color='purple', fontsize=8, ha='right')
            
            self.plot_hover_points.append({'x': x_curr, 'y': y_curr, 'label': f"S{num_plates} Top-L", 'type': 'stage'})
            self.plot_hover_points.append({'x': x_on_equilibrium, 'y': y_curr, 'label': f"S{num_plates} Top-R (on Eq.)\nx={x_on_equilibrium:.3f}", 'type': 'stage'})

            if self.feed_plate_number == "N/A" and x_on_equilibrium <= self.xi + self.CLS_TOLERANCE_EQ:
                 self.feed_plate_number = num_plates

            x_curr = x_on_equilibrium 
            
            if x_curr <= xb + self.CLS_TOLERANCE_EQ: 
                break 

            if x_curr > self.xi + self.CLS_TOLERANCE_EQ: 
                y_next_op = (R_val / (R_val + 1)) * x_curr + xd / (R_val + 1)
            else: 
                y_next_op = ss * (x_curr - xb) + xb
            
            if y_next_op < xb - self.CLS_TOLERANCE_EQ: 
                y_next_op = xb

            self.ax.plot([x_curr, x_curr], [y_curr, y_next_op], 'm-', lw=0.8)
            self.plot_hover_points.append({'x': x_curr, 'y': y_next_op, 'label': f"S{num_plates} Bot-R\ny={y_next_op:.3f}", 'type': 'stage'})
            
            y_curr = y_next_op
            if y_curr <= xb + self.CLS_TOLERANCE_EQ: 
                break

        if num_plates >= self.CLS_MAX_PLATES:
            if not first_run: messagebox.showwarning("Warning", f"Reached maximum number of plates ({self.CLS_MAX_PLATES}). Calculation might be stuck or require many stages.")

        self.current_num_plates = num_plates if num_plates < self.CLS_MAX_PLATES else f">= {self.CLS_MAX_PLATES}"
        self.num_plates_label.config(text=f"Number of plates: {self.current_num_plates}")
        self.feed_plate_label.config(text=f"Feed plate number: {self.feed_plate_number}")
        self.xi_label.config(text=f"Intersection xi: {self.xi:.4f}" if not np.isnan(self.xi) else "Intersection xi: --")
        self.yi_label.config(text=f"Intersection yi: {self.yi:.4f}" if not np.isnan(self.yi) else "Intersection yi: --")
        
        self.current_params_for_export = params 

        self._setup_plot_labels_and_legend()
        self.canvas.draw()

    def _setup_plot_labels_and_legend(self):
        self.ax.set_xlabel("X (mole fraction liquid)", fontsize=10)
        self.ax.set_ylabel("Y (mole fraction vapor)", fontsize=10)
        self.ax.set_xlim(0, 1)
        self.ax.set_ylim(0, 1)
        self.ax.legend(fontsize=8, loc='best')
        self.ax.grid(True, linestyle=':', alpha=0.7)
        self.ax.set_aspect('equal', adjustable='box')

    def export_plot(self):
        if not hasattr(self, 'current_params_for_export'):
            messagebox.showinfo("Export Plot", "Please calculate first to generate a plot.")
            return

        filepath = filedialog.asksaveasfilename(
            defaultextension=".png",
            filetypes=[("PNG files", "*.png"), ("JPEG files", "*.jpg"), ("PDF files", "*.pdf"), ("All files", "*.*")]
        )
        if not filepath: return

        fig_export, ax_export = plt.subplots(figsize=(8,7)) 
        params = self.current_params_for_export
        zf, q_val, xd, xb, R_val, alpha = params['zf'], params['q'], params['xd'], params['xb'], params['R'], params['alpha']

        x_eq_pts = np.linspace(0, 1, 200)
        y_eq_pts = self.equilib_y_from_x(x_eq_pts, alpha)
        ax_export.plot(x_eq_pts, y_eq_pts, 'r-', lw=1.5, label="Equilibrium Curve")
        ax_export.plot([0, 1], [0, 1], 'k-', lw=1, label="y=x")

        xi_exp, yi_exp = self.xi, self.yi 
        if np.isnan(xi_exp) or np.isnan(yi_exp): 
            if abs(q_val - 1.0) < self.CLS_TOLERANCE_EQ:
                xi_exp = zf
                yi_exp = (R_val / (R_val + 1)) * xi_exp + xd / (R_val + 1)
            else:
                m_q = q_val / (q_val - 1)
                c_q = zf * (1 - m_q)
                m_rol = R_val / (R_val + 1)
                c_rol = xd / (R_val + 1)
                if abs(m_q - m_rol) < self.CLS_TOLERANCE_DIV: 
                    messagebox.showerror("Export Error", "Cannot calculate intersection for export.")
                    plt.close(fig_export)
                    return
                xi_exp = (c_rol - c_q) / (m_q - m_rol)
                yi_exp = m_rol * xi_exp + c_rol
        
        ax_export.plot([xd, xi_exp], [xd, yi_exp], 'b-', lw=1, label="ROL")
        ss = (yi_exp - xb) / (xi_exp - xb) if abs(xi_exp-xb) > self.CLS_TOLERANCE_DIV else np.inf
        ax_export.plot([xb, xi_exp], [xb, yi_exp], 'b-', lw=1, label="SOL")
        ax_export.plot([zf, xi_exp], [zf, yi_exp], 'b--', lw=1, label="q-line")
        
        num_plates_exp = 0
        feed_plate_exp_local = "N/A" 

        x_curr, y_curr = xd, xd
        while x_curr > xb + self.CLS_TOLERANCE_EQ and num_plates_exp < self.CLS_MAX_PLATES:
            num_plates_exp += 1
            
            x_on_equilibrium_export = self.equilib_x_from_y(y_curr, alpha)
            
            ax_export.plot([x_curr, x_on_equilibrium_export], [y_curr, y_curr], 'm-', lw=0.8)
            ax_export.text(x_on_equilibrium_export - 0.01, y_curr + 0.02, str(num_plates_exp), color='purple', fontsize=8, ha='right')
            
            if feed_plate_exp_local == "N/A" and x_on_equilibrium_export <= xi_exp + self.CLS_TOLERANCE_EQ:
                 feed_plate_exp_local = num_plates_exp
            
            x_curr = x_on_equilibrium_export 
            
            if x_curr <= xb + self.CLS_TOLERANCE_EQ: 
                break 
            
            if x_curr > xi_exp + self.CLS_TOLERANCE_EQ: 
                y_next_op = (R_val / (R_val + 1)) * x_curr + xd / (R_val + 1)
            else: 
                y_next_op = ss * (x_curr - xb) + xb
            
            if y_next_op < xb - self.CLS_TOLERANCE_EQ: 
                y_next_op = xb

            ax_export.plot([x_curr, x_curr], [y_curr, y_next_op], 'm-', lw=0.8)
            y_curr = y_next_op
            
            if y_curr <= xb + self.CLS_TOLERANCE_EQ: 
                break
        
        num_plates_to_display = num_plates_exp if num_plates_exp < self.CLS_MAX_PLATES else f">= {self.CLS_MAX_PLATES}"

        ax_export.set_xlabel("X (mole fraction liquid)", fontsize=10)
        ax_export.set_ylabel("Y (mole fraction vapor)", fontsize=10)
        ax_export.set_xlim(0, 1); ax_export.set_ylim(0, 1)
        ax_export.legend(fontsize=8, loc='lower right', bbox_to_anchor=(1, 0.32))
        ax_export.grid(True, linestyle=':', alpha=0.7)
        ax_export.set_aspect('equal', adjustable='box')
        ax_export.set_title("McCabe-Thiele Diagram", fontsize=14)

        param_text = (
            f"zf = {params['zf']:.2f}\n"
            f"q = {params['q']:.2f}\n"
            f"xd = {params['xd']:.2f}\n"
            f"xb = {params['xb']:.2f}\n"
            f"R = {params['R']:.2f}\n"
            f"alpha = {params['alpha']:.2f}\n\n"
            f"Plates = {num_plates_to_display}\n"
            f"Feed Plate = {feed_plate_exp_local}\n"
            f"xi = {xi_exp:.3f}, yi = {yi_exp:.3f}"
        )
        fig_export.text(0.73, 0.05, param_text, 
             bbox=dict(boxstyle='round,pad=0.5', fc='lightyellow', alpha=0.8),
             fontsize=9, verticalalignment='bottom')

        fig_export.tight_layout(rect=[0, 0, 0.95, 1])
        
        try:
            fig_export.savefig(filepath)
            messagebox.showinfo("Export Plot", f"Plot saved to {filepath}")
        except Exception as e:
            messagebox.showerror("Export Error", f"Could not save plot: {e}")
        plt.close(fig_export)

    def _on_motion(self, event):
        if event.inaxes == self.ax:
            mouse_x_data, mouse_y_data = event.xdata, event.ydata
            if mouse_x_data is None or mouse_y_data is None: 
                if self.hover_text.get_visible():
                    self.hover_text.set_visible(False)
                    self.canvas.draw_idle()
                return

            found_point_info = None
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            min_dist_sq = ( (xlim[1]-xlim[0]) * self.CLS_HOVER_DIST_FACTOR )**2

            for point_info in self.plot_hover_points:
                if 'x' not in point_info or 'y' not in point_info:
                    continue
                dist_sq = (point_info['x'] - mouse_x_data)**2 + (point_info['y'] - mouse_y_data)**2
                if dist_sq < min_dist_sq:
                    min_dist_sq = dist_sq
                    found_point_info = point_info
            
            if found_point_info:
                self.hover_text.set_text(found_point_info['label'])
                offset_x = (xlim[1]-xlim[0]) * 0.015
                offset_y = (ylim[1]-ylim[0]) * 0.015
                pos_x = found_point_info['x'] + offset_x
                pos_y = found_point_info['y'] + offset_y
                
                try:
                    # Use a more robust way to estimate text box size in data coordinates
                    # This requires drawing to get extent, which can be slow if done every motion.
                    # For simplicity, we'll keep it this way or use an approximate fixed offset.
                    # A more advanced solution might cache text extents or use simpler logic.
                    bbox = self.hover_text.get_window_extent(renderer=self.canvas.get_renderer())
                    
                    # Transform bbox corners from display to data coordinates
                    inv = self.ax.transData.inverted()
                    bbox_data_bottom_left = inv.transform((bbox.x0, bbox.y0))
                    bbox_data_top_right = inv.transform((bbox.x1, bbox.y1))
                    
                    hover_width_data = abs(bbox_data_top_right[0] - bbox_data_bottom_left[0])
                    hover_height_data = abs(bbox_data_top_right[1] - bbox_data_bottom_left[1])
                    
                    if pos_x + hover_width_data > xlim[1]:
                        pos_x = found_point_info['x'] - offset_x - hover_width_data
                    if pos_y + hover_height_data > ylim[1]:
                        pos_y = found_point_info['y'] - offset_y - hover_height_data
                    if pos_x < xlim[0]: pos_x = xlim[0] # Avoid negative offset for left boundary
                    if pos_y < ylim[0]: pos_y = ylim[0] # Avoid negative offset for bottom boundary

                except Exception: 
                    pass # If renderer not ready or other issue, use default offset


                self.hover_text.set_position((pos_x, pos_y))
                if not self.hover_text.get_visible():
                    self.hover_text.set_visible(True)
                self.canvas.draw_idle()
            else:
                if self.hover_text.get_visible():
                    self.hover_text.set_visible(False)
                    self.canvas.draw_idle()
        else: 
            if self.hover_text.get_visible():
                self.hover_text.set_visible(False)
                self.canvas.draw_idle()

    def _save_parameters(self):
        params_to_save = {}
        for key, entry_widget in self.inputs.items():
            params_to_save[key] = entry_widget.get()
        
        filepath = filedialog.asksaveasfilename(
            defaultextension=".json",
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            title="Save Parameters"
        )
        if not filepath: return
        
        try:
            with open(filepath, 'w') as f:
                json.dump(params_to_save, f, indent=4)
            messagebox.showinfo("Save Parameters", f"Parameters saved to {filepath}")
        except Exception as e:
            messagebox.showerror("Save Error", f"Could not save parameters: {e}")

    def _load_parameters(self):
        filepath = filedialog.askopenfilename(
            filetypes=[("JSON files", "*.json"), ("All files", "*.*")],
            title="Load Parameters"
        )
        if not filepath: return
            
        try:
            with open(filepath, 'r') as f:
                loaded_params = json.load(f)
            
            for key, value in loaded_params.items():
                if key in self.inputs:
                    self.inputs[key].delete(0, tk.END)
                    self.inputs[key].insert(0, str(value))
            messagebox.showinfo("Load Parameters", f"Parameters loaded from {filepath}.\nPress Calculate to update plot.")
            # Optionally, trigger calculation automatically:
            # self.calculate_and_plot() 
        except Exception as e:
            messagebox.showerror("Load Error", f"Could not load parameters: {e}")
            
    def _reset_to_defaults(self):
        for key, entry_widget in self.inputs.items():
            entry_widget.delete(0, tk.END)
            entry_widget.insert(0, self.default_params[key])
        self.calculate_and_plot() 

    def show_about(self):
        about_text = """
This GUI demonstrates the McCabe-Thiele Graphical Method.
by Claudio A. Gelmi (2025) @ github.com/cagelmi

This is an educational tool for chemical engineering students.
It allows users to input distillation parameters and visualize the McCabe-Thiele diagram.

Reference:
McCabe, Smith and Harriott. Unit Operations of Chemical Engineering.

Hover over plot elements for more information.
        """
        messagebox.showinfo("About MCT", about_text)

if __name__ == '__main__':
    root = tk.Tk()
    app = McCabeThieleApp(root)
    root.mainloop()