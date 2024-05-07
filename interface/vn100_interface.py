from vnpy import *
#import Serial




if __name__ == "__main__":
	
	# declare sensor object
	s = VnSensor()
	print(type(s))
	
	# attempt connection: baud rate is bugged
	ez = EzAsyncData.connect('/dev/ttyUSB0', 115200)
	
	# alternate connection
	#s.connect('/dev/ttyUSB0', 9600)
	
	# info from sensor object
	print(s.is_connected)
	print(s.port)
	
	
	s.change_baudrate(9600)
	#print(s.read_serial_number())
	
	#print(s.port)
	print(s.read_model_number())
