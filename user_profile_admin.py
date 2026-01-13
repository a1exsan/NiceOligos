from nicegui import app, ui

class user_admin_profile():
    def __init__(self):

        with ui.column():
            with ui.row():
                ui.input(label="Пользователь", value=f'{app.storage.user.get("email")}')
                ui.input(label='Старый пароль:', password=True)
            with ui.row():
                ui.input(label='Новый пароль:', password=True)
                ui.input(label='Повторить новый пароль:', password=True)
            ui.button(text='Сохранить новый пароль')


