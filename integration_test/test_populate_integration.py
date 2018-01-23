import django.test
import populate


class TestPopulate(django.test.TestCase):

    def test_populate(self):
        populate.main('phageParser@example.com', 'integration_test/data', 5)
