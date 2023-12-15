
import turtle

# Initialize the graphics window
turtle.setup(width=800, height=600)
window = turtle.Screen()
window.title("Test")

# Set up the graphics environment
canvas = window.getcanvas()
canvas.config(width=800, height=600)

# Set up the turtle graphics
turtle.speed(0)
turtle.hideturtle()

# Call the Interface_Difr function
NameCatalogue = ''
Interface_Difr(NameCatalogue)

# Keep the graphics window open until it is closed by the user
turtle.mainloop()