;; We will just use the inherited init

FUNCTION postgres_adatc::write_status
  return,self.write_status
END 

PRO postgres_tsobj::tabledef, sqlfile

END


FUNCTION postgres_tsobj::input_file
  return,self->postgres_stuff_camcol::input_file()
END 

PRO postgres_tsobj::write_input

END 

PRO postgres_tsobj__define

  struct = {$
             postgres_tsobj, $
             write_status: 0, $
             INHERITS postgres_stuff_camcol $
           }

END 
