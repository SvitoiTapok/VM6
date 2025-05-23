import math
import tkinter as tk
from contextlib import nullcontext
from functools import reduce
from tkinter import ttk,filedialog, messagebox
import sympy as sp

import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import numpy as np
from functools import reduce
from numpy.lib.format import read_magic
from numpy.ma.core import zeros
from packaging.utils import InvalidName
from sympy.codegen.cfunctions import isnan
from sympy.physics.units import yards



def extra_Eiler(func, h,  x0, y0, b, accuracy, func_toch):
    out_data = [[x0], [y0]]
    x_i = x0
    y_i = y0
    while True:
        x_i1 = x_i + h
        if x_i1 > b+0.001:
            break
        y_tilda = y_i + h * func(x_i, y_i)
        y_i1 = y_i + h / 2 * (func(x_i, y_i) + func(x_i1, y_tilda))
        # print(y_i1, y_i, h, func(x_i, y_i), x_i, func(x_i1, y_tilda))
        x_i = x_i1
        y_i = y_i1
        out_data[1].append(y_i1)
        out_data[0].append(x_i1)
    out_data2 = [[x0], [y0]]
    h /= 2
    while True:
        print(h, "e")
        x_i = x0
        y_i = y0
        while True:
            x_i1 = x_i + h
            if x_i1 > b+0.001:
                break
            y_tilda = y_i + h * func(x_i, y_i)
            y_i1 = y_i + h / 2 * (func(x_i, y_i) + func(x_i1, y_tilda))
            # print(y_i1, y_i, h, func(x_i, y_i), x_i, func(x_i1, y_tilda))
            x_i = x_i1
            y_i = y_i1
            out_data2[1].append(y_i1)
            out_data2[0].append(x_i1)
        # print(h)
        # print(abs(out_data[1][-1] - out_data2[1][-1])/3)
        print(abs(out_data[1][-1] - out_data2[1][-1])/3)
        if abs(out_data[1][-1] - out_data2[1][-1])/3>accuracy:
            h /= 2
            out_data = out_data2
            out_data2 = [[x0], [y0]]
            continue
        return out_data2, h

def Runge(func, h,  x0, y0, b, accuracy, func_toch):
    out_data = [[x0], [y0]]
    x_i = x0
    y_i = y0
    while True:
        x_i1 = x_i + h
        if x_i1 > b+0.001:
            break
        k1 = h * func(x_i, y_i)
        k2 = h * func(x_i+h/2, y_i+k1/2)
        k3 = h * func(x_i+h/2, y_i+k2/2)
        k4 = h * func(x_i+h, y_i+k3)
        y_i1 = y_i + 1 / 6 * (k1+2*k2+2*k3+k4)
        x_i = x_i1
        y_i = y_i1
        out_data[1].append(y_i1)
        out_data[0].append(x_i1)
    out_data2 = [[x0], [y0]]
    h /= 2

    while True:
        print(h, "r")
        x_i = x0
        y_i = y0
        while True:
            x_i1 = x_i + h
            if x_i1 > b+0.001:
                break
            k1 = h * func(x_i, y_i)
            k2 = h * func(x_i + h / 2, y_i + k1 / 2)
            k3 = h * func(x_i + h / 2, y_i + k2 / 2)
            k4 = h * func(x_i + h, y_i + k3)
            y_i1 = y_i + 1 / 6 * (k1 + 2 * k2 + 2 * k3 + k4)
            x_i = x_i1
            y_i = y_i1
            out_data2[1].append(y_i1)
            out_data2[0].append(x_i1)
        if abs(out_data[1][-1] - out_data2[1][-1])/15>accuracy:
            h /= 2
            out_data = out_data2
            out_data2 = [[x0], [y0]]
            continue
        return out_data2, h


def Monnesy(func, h, x0, y0, b, accuracy, func_toch):


    while True:
        print(h, "m")
        data, h = Runge(func, h, x0, y0, x0 + 3 * h, accuracy, func_toch)
        y1, y2, y3, y4 = data[1][0], data[1][1], data[1][2], data[1][3]
        x1, x2, x3, x4 = data[0][0], data[0][1], data[0][2], data[0][3]
        out_data = [data[0][:4], data[1][:4]]
        while True:
            x_i = x4+h
            if x_i > b+0.001:
                break
            y_prog = y1 + 4*h/3*(2*func(x2, y2)-func(x3, y3)+2*func(x4, y4))
            y_kor = y3 + h/3*(func(x3, y3)+4*func(x4, y4)+func(x_i, y_prog))
            x1, x2, x3, x4 = x2, x3, x4, x_i
            y1, y2, y3, y4 = y2, y3, y4, y_kor
            out_data[0].append(x_i)
            out_data[1].append(y_kor)
        if max([abs(out_data[1][i]-func_toch(out_data[0][i], x0, y0)) for i in range(len(out_data[0]))]) > accuracy:
            h /= 2
            continue
        return out_data, h

def func1(x, y):
    return x*y
def real_func_1(x, x0, y0):
    C = y0/np.exp(x0**2/2)
    return C*np.exp(x**2/2)

def func2(x, y):
    return 2*y/x
def real_func_2(x, x0, y0):
    C = y0/x0**2
    return C*x**2

def func3(x, y):
    return -x*y**2-y
def real_func_3(x, x0, y0):
    C = (1/y0+x0+1)/np.exp(x0)
    return 1/(C*np.exp(x)-x-1)


# data, h = Monnesy(func3, 0.1, 1, 1, 3, 0.01, real_func_3)
# print(data)
#
# print(real_func_3(np.linspace(1, 1.5, int(np.ceil(0.5/h))), 1, 1))
# plt.scatter(data[0], data[1], color="red")
# plt.scatter(np.linspace(1, 3, int(np.ceil(0.5/h))), real_func_3(np.linspace(1, 3, int(np.ceil(0.5/h))), 1, 1), color="blue")
# plt.show()

# def get_a():
#     try:
#         a = float(x.replace(',', '.'))


def get_x0():
    try:
        return float(x0_entry.get().replace(',', '.'))
    except:
        messagebox.showerror("Ошибка!",
                             "Сообщение об ошибке: Ожидалось правильное значение x0")
    return None
def get_y0():
    try:
        return float(y0_entry.get().replace(',', '.'))
    except:
        messagebox.showerror("Ошибка!",
                             "Сообщение об ошибке: Ожидалось правильное значение y0")
    return None
def get_h():
    try:
        h = float(h_entry.get().replace(',', '.'))
        if h<=0:
            raise Exception
        return h
    except:
        messagebox.showerror("Ошибка!",
                             "Сообщение об ошибке: Ожидалось правильное значение h(число > 0)")
    return None
def get_xn():
    try:
        xn = float(int_label.get().replace(',', '.'))
    except:
        messagebox.showerror("Ошибка!",
                             "Сообщение об ошибке: Ожидалось правильное значение правой границы интервала")
        return None
    try:
        if float(x0_entry.get().replace(',', '.'))>xn:
            messagebox.showerror("Ошибка!",
                                 "Сообщение об ошибке: Ожидалось что правая граница будет больше левой")
            return None
        return xn
    except:
        return 1
def get_accuracy():
    try:
        accuracy = float(accuracy_entry.get().replace(',', '.'))
        if accuracy<=0:
            raise Exception
        return accuracy
    except:
        messagebox.showerror("Ошибка!",
                             "Сообщение об ошибке: Ожидалось правильное значение точности(число > 0")
    return None


# def read_data():
#     try:
#         data = input_area.get("1.0", tk.END)
#         data = list(set(data.strip().split('\n')))
#         d = [[],[]]
#         for i in data:
#             x, y = i.split()
#             d[0].append()
#             d[1].append(float(y.replace(',', '.')))
#         data = np.array(d)
#         if data.shape[1]<3:
#             delete_pre_res_func()
#             messagebox.showerror("Ошибка!","Сообщение об ошибке: Ожидалось хотя бы 3 разных точки для построения аппроксимирующей функции")
#             return 0
#         if len(set(d[0]))!=len(d[0]):
#             delete_pre_res_func()
#             messagebox.showerror("Ошибка!",
#                                  "Сообщение об ошибке: Ожидалось что будут переданы различные x")
#             return 0
#         sorted_data = list(zip(*sorted(zip(*data), key=lambda col: col[0])))
#         return sorted_data
#     except:
#         delete_pre_res_func()
#         messagebox.showerror("Ошибка!", "Сообщение об ошибке: Некорректные данные! Пожалуйста, введите данные в формате:\nx1 y1\nx2 y2")
#         return 0

def draw_func():
    global data_for_save
    delete_pre_res_func()
    x0, y0, h, xn, accuracy = get_x0(), get_y0(), get_h(), get_xn(), get_accuracy()
    if x0 is None or y0 is None or h is None or xn is None or accuracy is None:
        return 0
    ax.clear()
    draw_dots()
    fs=[]
    try:
        fs.append(["Усовершенствованный метод Эйлера", extra_Eiler(cur_func, h, x0, y0, xn, accuracy, cur_real_func), "red"])
        fs.append(["Метод Рунге-Кутта 4-го порядка", Runge(cur_func, h, x0, y0, xn, accuracy, cur_real_func), "yellow"])
        fs.append(["Метод Милна", Monnesy(cur_func, h, x0, y0, xn, accuracy, cur_real_func), "green"])
    except Exception as e:
        messagebox.showerror("Ошибка!",
                             f"Сообщение об ошибке: Вычислительная ошибка: {e}")
        return 0
    for i in range(len(fs)):
        data, h = fs[i][1]
        message = ""
        message += "Метод: " + fs[i][0] + '\n'
        message += "Исходные значения: функция, x0, y0, xn, точность: " + cur_func_name + str(x0) + ', ' + str(y0) + ', ' + str(xn) + ', ' + str(accuracy) + '\n'
        message += "Значения Х: " + str(data[0]) + '\n\n'
        message += "Рассчитанные значения Y: " + str(data[1]) + '\n\n'
        message += "Правильные значения Y: " + str([float(cur_real_func(data[0][i], x0, y0)) for i in range(len(data[0]))]) + '\n\n'
        message += "Полученное значение для h: " + str(h) + '\n\n'
        fs[i].append(message)
        ax.plot(data[0], data[1], color=fs[i][2], label=fs[i][0])
        #ax.scatter(data[0], data[1], color=fs[i][2], label=fs[i][0])
    fs = sorted(fs, key=lambda x: sum([abs(x[1][0][1][i]-cur_real_func(x[1][0][0][i], x0, y0)) for i in range(len(x[1][0]))]))
    data, h = fs[0][1]
    mes_for_root = ""
    mes_for_root += "Метод: " + fs[0][0] + '\n'
    mes_for_root += "Исходные значения: функция, x0, y0, xn, точность: " + cur_func_name + str(x0) + ', ' + str(y0) + ', ' + str(xn) + ', ' + str(accuracy) + '\n'
    mes_for_root += "Значения Х(первые 10): " + str(data[0][:10]) + '\n'
    mes_for_root += "Рассчитанные значения Y(первые 10): " + str(data[1][:10]) + '\n'
    mes_for_root += "Правильные значения Y(первые 10): " + str([float(cur_real_func(data[0][i], x0, y0)) for i in range(len(data[0]))][:10]) + '\n'
    mes_for_root += "Полученное значение для h: " + str(h) + '\n\n'
    root_label.config(text=mes_for_root)
    for i in fs:
        data_for_save += i[-1]

    # l = fs[0][1]
    # y_isk = 0
    # for k in range(len(l)):
    #     y_isk += l[len(l) - k - 1] * x_isk ** k
    # pol = ""
    # for k in range(len(l)):
    #     pol += f"({round(l[k], 3)})*x^{len(l) - k - 1}+"
    # pol = pol[:-5]
    # d_str = ""
    # for i in range(len(data[0])):
    #     d_str += str(data[0][i]) + " "
    # d_str+='\n'
    # for i in range(len(data[0])):
    #     d_str += str(data[1][i]) + " "
    # root_label.config(text=f"Входные параметры: \n{d_str}\nx={x_isk}, P(x) = {y_isk}, P = {pol}")
    # data_for_save = ""
    # for i in fs:
    #     l = i[1]
    #     y_isk = 0
    #     for k in range(len(l)):
    #         y_isk += l[len(l) - k - 1] * x_isk ** k
    #     for k in range(len(l)):
    #         pol += f"({round(l[k], 3)})*x^{len(l) - k - 1}+"
    #     pol = pol[:-5]
    #     data_for_save += f"Метод:{i[0]}\nВходные параметры: \n{d_str}\nx={x_isk}, P(x) = {y_isk}, P = {pol}\n"
    ax.legend(loc='upper right', fontsize=7, framealpha=1, shadow=True)
    ax.grid( color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
    ax.axhline(0, color='black', linewidth=0.5)
    canvas.draw()
    return 0

def draw_dots():
    global cur_real_func
    ax.clear()
    x0, y0, xn = get_x0(), get_y0(), get_xn()
    if x0 is None or y0 is None or xn is None:
        return 1
    x = np.linspace(x0, xn, 1000)
    y = cur_real_func(x, x0, y0)
    ax.plot(x, y,  color='blue', label='Реальная функция')
    ax.grid(color='gray', linestyle='--', linewidth=0.5, alpha=0.7)
    #ax.axhline(0, color='black', linewidth=0.5)
    canvas.draw()
    return 0



def delete_pre_res_func():
    global data_for_save
    Info_label.config(text="Информация о выполнении")
    root_label.config(text="")
    data_for_save=""


def show_frame(frame):
    frame.grid(row=0, column=0, sticky="nsew")

def save():
    if root_label.cget('text')=="":
        messagebox.showerror("Ошибка!",
                             f"Сообщение об ошибке: Перед сохранением рассчитайте значения функции, информацию о которых хотите сохранить")
        return 1
    filepath = filedialog.asksaveasfilename(
        defaultextension=".txt",
        filetypes=[("Текстовые файлы", "*.txt")]
    )
    if not filepath:
        return 1
    with open(filepath, "w", encoding="utf-8") as file:
        file.write(data_for_save)
    messagebox.showinfo("Успех", f"Файл сохранён: {filepath}")
def read_from_file():
    try:
        with open(file_name.get(), 'r') as file:
            data = file.read().strip()
            data = data.split()
            if len(data)!=5:
                raise ValueError
            x0_entry.delete(0, tk.END)
            x0_entry.insert(0, data[0])
            x0_entry.config(foreground='black')

            y0_entry.delete(0, tk.END)
            y0_entry.insert(0, data[1])
            y0_entry.config(foreground='black')

            int_label.delete(0, tk.END)
            int_label.insert(0, data[2])
            int_label.config(foreground='black')

            h_entry.delete(0, tk.END)
            h_entry.insert(0, data[3])
            h_entry.config(foreground='black')

            accuracy_entry.delete(0, tk.END)
            accuracy_entry.insert(0, data[4])
            accuracy_entry.config(foreground='black')

    except ValueError:
        messagebox.showerror("Ошибка!", f"Сообщение об ошибке: Некорректные данные")
    except FileNotFoundError:
        messagebox.showerror("Ошибка!", f"Сообщение об ошибке: Файл '{file_name.get()}' не найден")
    except Exception as e:
        messagebox.showerror("Ошибка!", f"Неизвестная ошибка: {e}")

def on_func_combobox_change(event):
    delete_pre_res_func()
    # Функция, которая вызывается при изменении значения в выпадающем списке
    global cur_func
    global cur_func_name
    global cur_real_func
    selected_value = func_combobox.get()  # Получаем выбранное значение
    if selected_value=="x*y":
        cur_func=func1
        cur_real_func = real_func_1
    elif selected_value=="2*y/x":
        cur_func=func2
        cur_real_func= real_func_2
    elif selected_value=="-x*y^2-y":
        cur_func=func3
        cur_real_func = real_func_3
    cur_func_name = selected_value
    messagebox.showinfo("Выбор функции", f"Вы выбрали: {selected_value}")  # Обновляем текст метки




def on_entry_click_x0(event):
    if x0_entry.get() == placeholder_text_x0:
        x0_entry.delete(0, tk.END)
        x0_entry.config(foreground='black')

def on_focusout_x0(event):
    if x0_entry.get() == '':
        x0_entry.insert(0, placeholder_text_x0)
        x0_entry.config(foreground='grey')

def on_entry_click_y0(event):
    if y0_entry.get() == placeholder_text_y0:
        y0_entry.delete(0, tk.END)
        y0_entry.config(foreground='black')

def on_focusout_y0(event):
    if y0_entry.get() == '':
        y0_entry.insert(0, placeholder_text_y0)
        y0_entry.config(foreground='grey')

def on_entry_click_int(event):
    if int_label.get() == placeholder_text_int:
        int_label.delete(0, tk.END)
        int_label.config(foreground='black')

def on_focusout_int(event):
    if int_label.get() == '':
        int_label.insert(0, placeholder_text_int)
        int_label.config(foreground='grey')

def on_entry_click_h(event):
    if h_entry.get() == placeholder_text_h:
        h_entry.delete(0, tk.END)
        h_entry.config(foreground='black')

def on_focusout_h(event):
    if h_entry.get() == '':
        h_entry.insert(0, placeholder_text_h)
        h_entry.config(foreground='grey')

def on_entry_click_accuracy(event):
    if accuracy_entry.get() == placeholder_text_accuracy:
        accuracy_entry.delete(0, tk.END)
        accuracy_entry.config(foreground='black')

def on_focusout_accuracy(event):
    if accuracy_entry.get() == '':
        accuracy_entry.insert(0, placeholder_text_accuracy)
        accuracy_entry.config(foreground='grey')

placeholder_text_x0 = "X0"
placeholder_text_y0 = "Y0"
placeholder_text_int = "Правая граница интервала x_n"
placeholder_text_h = "h"
placeholder_text_accuracy = "Точность"


cur_func=func1
cur_real_func = real_func_1
cur_func_name="sin(exp(x) + x)"
# cur_meth_name="Метод левых прямоугольников"

data_for_save=""

root = tk.Tk()
root.title("Динамическое построение графиков")
root.grid_rowconfigure(0, weight=1)
root.grid_columnconfigure(0, weight=1)

int_frame = ttk.Frame(root)
int_frame.grid_rowconfigure(0, weight=1)
int_frame.grid_rowconfigure(1, weight=1)
int_frame.grid_rowconfigure(2, weight=1)
int_frame.grid_rowconfigure(3, weight=1)
int_frame.grid_rowconfigure(4, weight=1)
int_frame.grid_rowconfigure(5, weight=1)
int_frame.grid_rowconfigure(6, weight=1)
int_frame.grid_rowconfigure(7, weight=1)
int_frame.grid_rowconfigure(8, weight=1)
int_frame.grid_rowconfigure(9, weight=1)
int_frame.grid_rowconfigure(10, weight=1)
int_frame.grid_rowconfigure(11, weight=1)
int_frame.grid_rowconfigure(12, weight=1)
int_frame.grid_columnconfigure(0, weight=1)
int_frame.grid_columnconfigure(1, weight=1)

# # Поле для ввода погрешности
# accuracy_label = ttk.Label(int_frame, text="Точность:")
# accuracy_label.grid(row=0,column=0, padx=10, pady=10,sticky="ew")
# accuracy = ttk.Entry(int_frame)
# accuracy.insert(0, "0.01")
# accuracy.grid(row=0,column=1, padx=10, pady=10,sticky="ew")
#
# # Поле для ввода отрезка
# left_gran_label = ttk.Label(int_frame, text="левая граница:")
# left_gran_label.grid(row=1,column=0, padx=10, pady=10,sticky="ew")
# left_gran = ttk.Entry(int_frame)
# left_gran.insert(0, "0.0")
# left_gran.grid(row=1,column=1, padx=10, pady=10,sticky="ew")
#
#
# right_gran_label = ttk.Label(int_frame, text="правая граница:")
# right_gran_label.grid(row=2,column=0, padx=10, pady=10,sticky="ew")
# right_gran = ttk.Entry(int_frame)
# right_gran.insert(0, "1.0")
# right_gran.grid(row=2,column=1, padx=10, pady=10,sticky="ew")
row = 1


x0_entry = ttk.Entry(int_frame, width=30)
x0_entry.grid(row=row, column=0, pady=10)

x0_entry.insert(0, placeholder_text_x0)
x0_entry.config(foreground='grey')
x0_entry.bind('<FocusIn>', on_entry_click_x0)
x0_entry.bind('<FocusOut>', on_focusout_x0)

y0_entry = ttk.Entry(int_frame, width=30)
y0_entry.grid(row=row, column=1, pady=10)

y0_entry.insert(0, placeholder_text_y0)
y0_entry.config(foreground='grey')
y0_entry.bind('<FocusIn>', on_entry_click_y0)
y0_entry.bind('<FocusOut>', on_focusout_y0)

row+=1

int_label = ttk.Entry(int_frame, width=30)
int_label.grid(row=row, column=0, pady=10)

int_label.insert(0, placeholder_text_int)
int_label.config(foreground='grey')
int_label.bind('<FocusIn>', on_entry_click_int)
int_label.bind('<FocusOut>', on_focusout_int)

h_entry = ttk.Entry(int_frame, width=30)
h_entry.grid(row=row, column=1, pady=10)

h_entry.insert(0, placeholder_text_h)
h_entry.config(foreground='grey')
h_entry.bind('<FocusIn>', on_entry_click_h)
h_entry.bind('<FocusOut>', on_focusout_h)

row+=1


functions = ["x*y", "2*y/x", "-x*y^2-y"]
func_combobox = ttk.Combobox(int_frame, values=functions)
func_combobox.set("x*y")
func_combobox.bind("<<ComboboxSelected>>", on_func_combobox_change)

func_combobox.grid(row=row,column=0, pady = 10)

accuracy_entry = ttk.Entry(int_frame, width=30)
accuracy_entry.grid(row=row, column=1, pady=10)

accuracy_entry.insert(0, placeholder_text_accuracy)
accuracy_entry.config(foreground='grey')
accuracy_entry.bind('<FocusIn>', on_entry_click_accuracy)
accuracy_entry.bind('<FocusOut>', on_focusout_accuracy)

row+=1

fig, ax = plt.subplots()
canvas = FigureCanvasTkAgg(fig, master=int_frame)
draw_button = ttk.Button(int_frame, text="Нарисовать точную функцию", command=draw_dots)
draw_button.grid(row=row, column=0, padx=10, pady=10,sticky="ew")
draw_button = ttk.Button(int_frame, text="Рассчитать приближенное решение", command=draw_func)
draw_button.grid(row=row, column=1, padx=10, pady=10,sticky="ew")
row+=1

save_button = ttk.Button(int_frame, text="Сохранить результат", command=save)
save_button.grid(row=row, columnspan=2, padx=10, pady=10,sticky="ew")
row+=1
canvas.get_tk_widget().grid(row=row, column=0, columnspan=2)
row+=1

Info_label = ttk.Label(int_frame, text="Информация о выполнении:")
Info_label.grid(row=row, column=0, columnspan=2)
row+=1


root_label = ttk.Label(int_frame, text="")
root_label.grid(row=row, column=0, columnspan=2)
row+=1


file_name_label = ttk.Label(int_frame, text="Название файла:")
file_name_label.grid(row=row, column=0, padx=10, pady=10, sticky="ew")

file_name = ttk.Entry(int_frame)
file_name.insert(0, "a.txt")
file_name.grid(row=row, column=1, padx=10, pady=10, sticky="ew")
row+=1

load_f_button = ttk.Button(int_frame, text="Загрузить из файла", command=read_from_file)
load_f_button.grid(row=row, column=0, columnspan=2)
row+=1

show_frame(int_frame)
# Запуск основного цикла
root.mainloop()