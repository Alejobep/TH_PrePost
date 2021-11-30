
from bokeh.models.glyphs import Line
from bokeh.models.selections import IntersectRenderers
from bokeh.palettes import Viridis256
from bokeh.plotting import figure
from bokeh.transform import  factor_cmap
from bokeh.models import LinearColorMapper
import numpy as np


def set_color(df, col, palette = Viridis256):
    """
    Sets the color settings tu use in Cross section plot
    datramapper bases on parameter

    """
    if df[col].apply(isinstance,args = [float]).all():
        mapper = LinearColorMapper(palette=Viridis256, 
                                   low=df[col].min(), 
                                   high=df[col].max())

        color = {'field': col, 'transform': mapper}

    else:

        categories = df[col].drop_duplicates().values

        if len(categories)==1:
            color_list = palette[0:1]

        else:
            color_list = palette[::(len(palette)-1)//(len(categories)-1)]

        color = factor_cmap(col,
                            palette=color_list, 
                            factors=categories)
    return color


def set_str_tooltips(df, *col):
    
    tt_text = [('cell','@ElName')]

    col_unit = {'P':'bar',
                'Pressure':'bar',
                'T':'degC',
                'Temperature':'degC',
                'S_hyd':'frac.',
                'S_Hydrate':'frac.',
                'S_gas':'frac.',
                'S_aqu':'frac.',
                'S_aqueous':'frac.',
                'S_icd':'frac.',
                'S_ice':'frac.',
                'X_inh':'frac.',
                'X_Inhibitor':'frac.',
                'C_CH4inGas':'Kg/m3',
                'C_CH4inAqu':'Kg/m3',
                'Dens_Gas':'Kg/m3',
                'Dens_Aqu':'Kg/m3',
                'Dens_Hydr':'Kg/m3',
                'Visc_Gas':'Kg/m*s',
                'Visc_Aqu':'Kg/m*s',
                'k_rg':'frac.',
                'Krel_Gas':'frac.',
                'k_rw':'frac.',
                'Krel_Aqu':'frac.',
                'k_adj_F':'frac.',
                'perm_abs':'mD',
                'porosity':'frac.',
                'P_CH4':'bar',
                'P_EqHydr':'bar',
                'P_SatWat':'bar',
                'P_cap':'bar',
                'P_aqu':'bar'}
    
    for c in col:
        col_name = c

        if df[c].apply(isinstance,args = [float]).all():

            col_value = '@'+c+' '+col_unit[c]

            tt_text.append((col_name,col_value))

        else:

            col_value = '@'+c


            tt_text.append((col_name,col_value))
    return tt_text
        
def create_figure1(df,col,time_steps,data_source, plane = 'xy'):

    #Define ColumnDataSource
    source_2D=data_source
    
    #Set X, Y ranges
    x=plane[0]
    y=plane[1]

    first_row=df[df[y]==df[y].max()].loc[time_steps[0]]
    last_row= df[df[y]==df[y].min()].loc[time_steps[0]]

    first_col= df[df[x]==df[x].min()].loc[time_steps[0]]
    last_col=df[df[x]==df[x].max()].loc[time_steps[0]]

    Xrange=(first_col[x].min()-first_col['d'+x].max()/2,last_col[x].max()+last_col['d'+x].max()/2)
    Yrange=(last_row[y].min()-last_row['d'+y].max()/2,first_row[y].max()+first_row['d'+y].max()/2)

    #Define
    
    plot_options = dict(toolbar_location=None,
                        y_range=Yrange,
                        x_range=Xrange)    


    figure1 = figure(width=500, x_axis_location='above', plot_height=350,
                     **plot_options)

    label_X = x+' distance [m]'
    label_Y = y+' distance [m]'


    figure1.xaxis.axis_label = label_X

    figure1.yaxis.axis_label = label_Y
    
    figure1.rect(x=x,y=y,width='d'+x,height='d'+y,source=source_2D,
         fill_color=set_color(df,col), line_color=None, line_width=0)  

    return figure1

def create_figure2(df, time_steps, data_source, plane ='xy' ):

    #Define ColumnDataSource
    source_2D=data_source

    #Set X, Y ranges
    x=plane[0]
    y=plane[1]

    first_row=df[df[y]==df[y].max()].loc[time_steps[0]]
    last_row= df[df[y]==df[y].min()].loc[time_steps[0]]

    first_col= df[df[x]==df[x].min()].loc[time_steps[0]]
    last_col=df[df[x]==df[x].max()].loc[time_steps[0]]

    Xrange=(first_col[x].min()-first_col['d'+x].max()/2,last_col[x].max()+last_col['d'+x].max()/2)
    Yrange=(last_row[y].min()-last_row['d'+y].max()/2,first_row[y].max()+first_row['d'+y].max()/2)


    plot_options = dict(toolbar_location=None,
                        y_range=Yrange,
                        x_range=Xrange)    


    figure2 = figure(width=500, x_axis_location='below', plot_height=350,
                    #  tooltips = TOOLTIPS,**plot_options)
                    **plot_options)

    label_X = x+' distance [m]'

    if y == 'z':
        label_Y = 'depth [m]'

    else:
        label_Y = y+' distance [m]'



    figure2.xaxis.axis_label = label_X

    figure2.yaxis.axis_label = label_Y

    #Draw brine
    figure2.rect(x=x,y=y,width='d'+x,height='d'+y,source=source_2D,
        fill_color="#AB957C",fill_alpha={'field':'S_aqu'}, line_color="#AB957C", line_width = 0.5, line_alpha={'field':'S_aqu'}) 

    #Draw gas
    figure2.rect(x=x,y=y,width='d'+x,height='d'+y,source=source_2D,
        fill_color="firebrick",fill_alpha={'field':'S_gas'},line_color="firebrick", line_width = 0.5, line_alpha={'field':'S_gas'}) 



    
    return figure2

def create_figure1_1D(df,col,time_steps,data_source, axis = 'y'):
    #Define ColumnDataSource
    source_1D=data_source

    #Set X, Y ranges
    x=col
    y=axis

    first_row=df[df[y]==df[y].max()].loc[time_steps[0]]
    last_row= df[df[y]==df[y].min()].loc[time_steps[0]]

    if df[col].apply(isinstance,args = [float]).all():
        minv = df[col].min()
        maxv = df[col].max()

        if minv/maxv > 0.97:
            delta = (10**np.floor(np.log10(np.abs(minv))))
            minv -= delta
            maxv += delta

        Xrange= (minv, maxv)

    else:
        Xrange = (0,1)

    Yrange=(last_row[y].min()-last_row['d'+y].max()/2,first_row[y].max()+first_row['d'+y].max()/2)

    #Define

    plot_options = dict(toolbar_location=None,
                        y_range=Yrange, x_range=Xrange)    


    figure1_1D = figure(width=170, x_axis_location='above', y_axis_location='right', plot_height=350,
                    **plot_options)

    label_X = x
    label_Y = axis+' distance [m]'


    figure1_1D.xaxis.axis_label = label_X

    figure1_1D.yaxis.axis_label = label_Y

    figure1_1D.line(x=x, y=y,  line_width=2, color="indigo",source=source_1D)    

    return figure1_1D

def create_figure2_1D(df,time_steps,data_source, axis = 'y'):
    #Define ColumnDataSource
    source_1D=data_source

    #Set X, Y ranges
    y=axis

    first_row=df[df[y]==df[y].max()].loc[time_steps[0]]
    last_row= df[df[y]==df[y].min()].loc[time_steps[0]]


    Yrange=(last_row[y].min()-last_row['d'+y].max()/2,first_row[y].max()+first_row['d'+y].max()/2)

    #Define

    plot_options = dict(toolbar_location=None,
                        y_range=Yrange)    


    figure2_1D = figure(width=170, x_axis_location='below', y_axis_location='right', plot_height=350,
                    **plot_options)

    label_X = 'saturation [frac.]'
    label_Y = axis+' distance [m]'


    figure2_1D.xaxis.axis_label = label_X

    figure2_1D.yaxis.axis_label = label_Y
    #Draw gas

    figure2_1D.line(x='S_aqu', y=y,  line_width=2, color="steelblue",source=source_1D)    
    figure2_1D.line(x='S_gas', y=y,  line_width=2, color="firebrick",source=source_1D)    
    figure2_1D.line(x='S_hyd', y=y,  line_width=2, color="black",source=source_1D)    

    return figure2_1D

